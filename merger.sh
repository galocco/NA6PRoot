
DIR=tof_40ag_ev_abs

rm $DIR/*root

hadd -f HitsMuonSpecModular.root $DIR/job_*/HitsMuonSpecModular.root
hadd -f HitsVerTel.root $DIR/job_*/HitsVerTel.root
hadd -f MCKine.root $DIR/job_*/MCKine.root

hadd -f ClustersMuonSpec.root $DIR/job_*/ClustersMuonSpec.root
hadd -f ClustersTOF.root $DIR/job_*/ClustersTOF.root
hadd -f ClustersVerTel.root $DIR/job_*/ClustersVerTel.root
hadd -f HitsMuonSpecModular.root $DIR/job_*/HitsMuonSpecModular.root
hadd -f HitsTOF.root $DIR/job_*/HitsTOF.root
hadd -f HitsVerTel.root $DIR/job_*/HitsVerTel.root
hadd -f MCKine.root $DIR/job_*/MCKine.root
hadd -f TracksTOFMatching.root $DIR/job_*/TracksTOFMatching.root
hadd -f TracksVerTel.root $DIR/job_*/TracksVerTel.root
hadd -f VerticesVerTel.root $DIR/job_*/VerticesVerTel.root

cp $DIR/job_0/na6pLayout.ini $DIR/
cp $DIR/job_0/na6pRecoParam.ini $DIR/
cp $DIR/job_0/geometry.root $DIR/

cp ClustersMuonSpec.root $DIR/
cp ClustersTOF.root $DIR/
cp ClustersVerTel.root $DIR/
cp HitsMuonSpecModular.root $DIR/
cp HitsTOF.root $DIR/
cp HitsVerTel.root $DIR/
cp MCKine.root $DIR/
cp TracksTOFMatching.root $DIR/
cp TracksVerTel.root $DIR/
cp VerticesVerTel.root $DIR/