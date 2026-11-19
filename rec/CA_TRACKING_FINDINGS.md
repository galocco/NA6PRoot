# CA tracking: propagation reuse and cluster ordering

## Current implementation

The implementation is in [NA6PTrackerCA.cxx](src/NA6PTrackerCA.cxx), with
declarations in [NA6PTrackerCA.h](include/NA6PTrackerCA.h).

- `computeTrackToClusterChi2` always propagates the supplied cell track in place
  to the cluster's actual Z and then evaluates its predicted chi-square.
- The `caPropagateCellTrackForChi2` configuration option, tracker flag, and setter
  have been removed. The user reported an improvement with the enabled behavior;
  this document does not assign a measured speedup or efficiency change to it.
- Both the pairwise chi-square cache in `prolongSeed` and the road extension cache
  in `findRoads` have been removed.
- The propagation savings report and its counters have been removed.
- Road extension defaults to inward and then outward. Set
  `reco.caDoInwardSeedExtension=false` to skip inward seed extension for both VT
  and MS and extend the original seeds outward only. This option does not control
  inward refitting or propagation to the primary vertex.
- Cluster sorting has been restored to `sortClustersByLayerAndY`: group by layer
  ID in linear time, then sort each layer by Y. Tracklet range searches operate
  directly on the cluster array. The separate Y index, Z range checks, and Z
  sorting have been removed. Fixed-size arrays for layer counts and Z sums remain.

## Z ordering experiment (reverted)

The following records the explored implementation and its findings. It is no
longer active; the current implementation uses layer-plus-Y sorting as described
above.

`sortClustersByZAndY` uses the relationship between layer ID and Z to avoid a
global cluster sort:

1. Count clusters and accumulate their Z sum, minimum, and maximum per layer.
2. Compute layer offsets and group clusters by layer ID in linear time.
3. Sort clusters within each layer by increasing Z, then increasing Y when Z is
   equal.
4. Supply a Y-ordered index for the tracklet range search.

There is a shortcut for layers whose clusters all have exactly the same Z:
sorting by Y alone already gives the required `(Z, Y)` order. Their Y lookup
indices remain in identity order, so there is no second sort. Empty and
single-cluster layers also need no sorting.

Only layers with varying Z need both the `(Z, Y)` cluster sort and a separate
sort of their Y lookup indices in `mClustersByY`. Tracklets and cells retain
their existing index-based sorting and matching logic.

The layer counts and Z statistics use fixed-size arrays for the existing
25-layer limit. There is no sort of layer means and no global cluster sort.

### Geometry assumption

For the grouped result to equal a global `(Z, Y)` sort, populated layers must
have disjoint Z ranges increasing with layer ID:

```text
maxZ(previous populated layer) < minZ(next populated layer)
```

The code checks this condition and reports a fatal error if it is violated.
Different layers having no identical Z values is not, by itself, sufficient:
their ranges must also not interleave. Clusters within one layer may have
different Z values, provided the condition above holds.

`layerID` is not part of the cluster comparison. It is used for efficient
grouping and for the `firstIndex`/`lastIndex` ranges required by CA.

## Why the separate Y index is necessary

`computeLayerTracklets` restricts candidates with `std::partition_point` on a
Y interval. That operation requires Y-sorted input. A layer sorted by `(Z, Y)`
is not necessarily sorted by Y when Z varies.

The search therefore operates on `mClustersByY` and retrieves each cluster using
the resulting index. This preserves the binary range lookup without changing
the cluster array's Z ordering.

## Performance implications

The original layer-plus-Y version grouped clusters in linear time and sorted
each layer by Y. The intermediate global `(Z, Y)` implementation added a global
sort and then sorted the per-layer Y indices as well.

The final version removes that global sort. For constant-Z layers, it performs
one local Y sort, as the original version did. Staggered layers incur an
additional local index sort. Grouping, statistics, and index initialization
remain linear in the number of clusters; the Y index stores one integer per
cluster, and grouping uses a temporary cluster vector.

The optimized experiment reduced sorting work relative to the intermediate
version. Its effect on total reconstruction time was not benchmarked. The
subsequent restoration of layer-plus-Y removes the separate Y index entirely
and uses a single local Y sort for every layer, including staggered layers.

## Limits of Z ordering for propagation reuse

Sorting clusters by Z does not guarantee that every cell track is propagated
monotonically. `findCellsNeighbours` and `prolongSeed` reuse mutable cell fits,
and different road candidates or seeds can revisit the same fit.

For example, one scan can request Z values `10, 11`, and a later scan using the
same cell can request `10, 11` again. The transition between scans still includes
`11 -> 10`, even though each scan is ordered.

No additional propagation-state cache or batching of chi-square requests was
introduced. A reduction in direction reversals has not been measured.

### Remaining staggered-geometry consideration

`mLayersZ` still stores the mean cluster Z per layer for estimating the tracklet
Y search window. It is not used to sort clusters. Exact candidate checks and
chi-square propagation use each cluster's actual Z.

For a sufficiently large Z spread within a layer, a window based on mean Z may
exclude a valid candidate before the exact check. Computing a conservative
window over the layer's Z range remains a separate, unimplemented change.

## Earlier propagation counters

The old `same-Z` and `road reuse` values described separate mechanisms:

- Same-Z savings counted checks that could directly evaluate chi-square at the
  track's current Z when cell-state reuse was enabled.
- Road reuse counted cache hits and skipped extension checks in the copy-based
  path. Those caches were bypassed when cell-state reuse was enabled.

This explained why `road reuse` could be positive with cell-track reuse disabled
and zero with it enabled. Those caches and reporting counters are now removed.

## Optional outward-only seed extension

The default `caDoInwardSeedExtension=true` preserves the existing search. Disable
it in the `[reco]` section of a reconstruction INI file:

```ini
[reco]
caDoInwardSeedExtension=false
```

Alternatively, append `reco.caDoInwardSeedExtension=false` to the command's
`--configKeyValues` argument (separate it from other settings with a semicolon).

For five layers with cells `(0,1,2)`, `(1,2,3)`, and `(2,3,4)`, the earliest
accepted neighbouring-cell pair seeds layers `0..3`; outward extension can add
layer 4. Inward extension is redundant for that successful path.

Equal reconstruction efficiency for outward-only extension is not established.
Later accepted seeds can currently recover earlier clusters through their own
extension checks. The option enables a separate comparison; equal efficiency has
not been measured.

## Verification

- `cmake --build build --target recLib -j5` completed successfully after restoring
  layer-plus-Y sorting.
- `git diff --check` passed.
- The inspected build has no registered CTest tests.
- No timing benchmark, propagation reversal measurement, or reconstruction
  efficiency comparison was run for the final ordering change.
- The optional inward seed extension change was built and installed locally.
  One event was reconstructed for both VT and MS with the option set to `true`
  and then `false`, using separate temporary output directories. Both runs
  completed successfully and produced one track per detector. This smoke check
  does not establish equal reconstruction efficiency on a larger sample.
