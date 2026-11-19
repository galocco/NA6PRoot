import configparser
import optuna
import subprocess
import ROOT
import logging
from pathlib import Path
from typing import Dict, Optional
import sys
import argparse

# Configuration
INI_FILE = "na6pMSCA.ini"
OUTPUT_ROOT_FILE = "MSTrackFinderCAoutput.root"
TRIAL_CONFIG_FILE = "trial_config.ini"
OPTIMIZED_CONFIG_FILE = "optimized_config.ini"
DB_PATH = "sqlite:///optuna_ms_ca.db"
STUDY_NAME = "ms_ca_optimization"

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('optimization.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class CaseSensitiveConfigParser(configparser.ConfigParser):
    """ConfigParser that preserves the case of keys"""
    def optionxform(self, optionstr):
        return optionstr


def validate_files() -> bool:
    """Check that required files exist"""
    if not Path(INI_FILE).exists():
        logger.error(f"INI file not found: {INI_FILE}")
        return False
    return True


def get_first_trial_from_ini() -> Dict:
    """
    Read current parameter values from INI file to use as first trial.
    """
    config = CaseSensitiveConfigParser()
    config.read(INI_FILE)
    
    first_trial = {}
    
    try:
        # Read number of iterations
        nIterationsTrackerCA = int(config["reco"]["nIterationsTrackerCA"])
        first_trial["nIterationsTrackerCA"] = nIterationsTrackerCA
        
        # Read all CA parameters for each iteration
        params_to_read = [
            "maxDeltaThetaTrackletsCA",
            "maxDeltaPhiTrackletsCA",
            "maxDeltaTanLCellsCA",
            "maxDeltaPhiCellsCA",
            "maxDeltaPxPzCellsCA",
            "maxDeltaPyPzCellsCA",
            "maxChi2TrClCellsCA",
            "maxChi2ndfCellsCA",
            "maxChi2ndfTracksCA",
            "minNClusTracksCA"
        ]
        
        for i in range(nIterationsTrackerCA):
            for param in params_to_read:
                key = f"{param}[{i}]"
                if key in config["reco"]:
                    value = config["reco"][key]
                    # Convert to appropriate type
                    if param.startswith("minNClus"):
                        first_trial[f"{param}{i}"] = int(value)
                    else:
                        first_trial[f"{param}{i}"] = float(value)
                else:
                    logger.warning(f"Parameter {key} not found in INI file")
        
        logger.info(f"Read {len(first_trial)} parameters from INI file")
        return first_trial
        
    except Exception as e:
        logger.error(f"Error reading INI file: {e}")
        raise


def objective(trial: optuna.Trial, n_events: int = 1000) -> float:
    """
    Objective function for Optuna.
    Suggests new values for parameters and returns a metric to optimize.
    """
    try:
        # Load config with case-sensitive parser
        config = CaseSensitiveConfigParser()
        config.read(INI_FILE)
        
        # Optimize number of iterations
        nIterationsTrackerCA = trial.suggest_int("nIterationsTrackerCA", 1, 4)
        config["reco"]["nIterationsTrackerCA"] = str(nIterationsTrackerCA)

        # Update config with new trial parameters
        for i in range(nIterationsTrackerCA):
            # Suggest parameters with appropriate ranges
            params = {
                f"maxDeltaThetaTrackletsCA{i}": trial.suggest_float(
                    f"maxDeltaThetaTrackletsCA{i}", 0.02, 1, log=True
                ),
                f"maxDeltaPhiTrackletsCA{i}": trial.suggest_float(
                    f"maxDeltaPhiTrackletsCA{i}", 0.05, 1, log=True
                ),
                f"maxDeltaTanLCellsCA{i}": trial.suggest_float(
                    f"maxDeltaTanLCellsCA{i}", 2, 20
                ),
                f"maxDeltaPhiCellsCA{i}": trial.suggest_float(
                    f"maxDeltaPhiCellsCA{i}", 0.1, 2.4
                ),
                f"maxDeltaPxPzCellsCA{i}": trial.suggest_float(
                    f"maxDeltaPxPzCellsCA{i}", 0.01, 0.1
                ),
                f"maxDeltaPyPzCellsCA{i}": trial.suggest_float(
                    f"maxDeltaPyPzCellsCA{i}", 0.01, 0.1
                ),
                f"maxChi2TrClCellsCA{i}": trial.suggest_float(
                    f"maxChi2TrClCellsCA{i}", 1, 20
                ),
                f"maxChi2ndfCellsCA{i}": trial.suggest_float(
                    f"maxChi2ndfCellsCA{i}", 1, 20
                ),
                f"maxChi2ndfTracksCA{i}": trial.suggest_float(
                    f"maxChi2ndfTracksCA{i}", 1, 20
                ),
            }
            
            # Write to config
            for param_key, param_value in params.items():
                # Extract base parameter name
                param_name = param_key.rstrip('0123456789')
                config["reco"][f"{param_name}[{i}]"] = str(param_value)
            
            config["reco"][f"minNClusTracksCA[{i}]"] = str(6)
        
        # Save updated config
        with open(TRIAL_CONFIG_FILE, "w") as f:
            config.write(f, space_around_delimiters=False)
        
        logger.info(f"Trial {trial.number}: Running ROOT macro...")
        
        # Run ROOT macro with trial config
        # NB: reduce number of events during tuning to avoid long hangs
        result = subprocess.run(
            ["root", "-l", "-b", "-q",
             f'runMSTrackFinderCA.C+(0, {n_events}, "../jpsi50gev", 6, "{TRIAL_CONFIG_FILE}")'],
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        if result.returncode != 0:
            logger.error(f"ROOT error (trial {trial.number}): {result.stderr}")
            return -1e6
        
        # Compute metric
        metric = compute_metric(trial.number)
        logger.info(f"Trial {trial.number}: Metric = {metric:.4f}")
        
        return metric
        
    except subprocess.TimeoutExpired:
        logger.error(f"Trial {trial.number}: Timeout during ROOT execution")
        return -1e6
    except Exception as e:
        logger.error(f"Trial {trial.number}: Unexpected error: {e}")
        return -1e6


def compute_metric(trial_number: Optional[int] = None) -> float:
    """
    Compute the metric from the reconstruction output.
    """
    try:
        if not Path(OUTPUT_ROOT_FILE).exists():
            logger.error(f"Output file not found: {OUTPUT_ROOT_FILE}")
            return -1e6
        
        file = ROOT.TFile.Open(OUTPUT_ROOT_FILE)
        if not file or file.IsZombie():
            logger.error("Unable to open ROOT file")
            return -1e6
        
        # Read metrics
        eff = file.Get("efficiency")
        purity = file.Get("purity")
        time = file.Get("timePerEvent")
        
        if not all([eff, purity, time]):
            logger.error("One or more metrics not found in ROOT file")
            file.Close()
            return -1e6
        
        eff_val = eff.GetVal()
        purity_val = purity.GetVal()
        time_val = time.GetVal()
        
        # Metric formula (adjustable as needed)
        # Favors efficiency, considers purity, penalizes time
        metric = eff_val + purity_val / 7.0# - time_val / 0.001
        
        log_msg = (f"{'Trial ' + str(trial_number) + ': ' if trial_number is not None else ''}"
                   f"Efficiency={eff_val:.4f}, Purity={purity_val:.4f}, "
                   f"Time={time_val:.4f}s, Metric={metric:.4f}")
        logger.info(log_msg)
        
        file.Close()
        return metric
        
    except Exception as e:
        logger.error(f"Error computing metric: {e}")
        return -1e6


def save_best_config(study: optuna.Study) -> None:
    """Save the best configuration to an INI file"""
    try:
        config = CaseSensitiveConfigParser()
        config.read(INI_FILE)
        
        trial = study.best_trial
        nIter = trial.params["nIterationsTrackerCA"]
        config["reco"]["nIterationsTrackerCA"] = str(nIter)
        
        param_names = [
            "maxDeltaThetaTrackletsCA", "maxDeltaPhiTrackletsCA", 
            "maxDeltaTanLCellsCA", "maxDeltaPhiCellsCA", 
            "maxDeltaPxPzCellsCA", "maxDeltaPyPzCellsCA",
            "maxChi2TrClCellsCA", "maxChi2ndfCellsCA", 
            "maxChi2ndfTracksCA", "minNClusTracksCA"
        ]
        
        for i in range(nIter):
            for param in param_names:
                key = f"{param}{i}"
                if key in trial.params:
                    config["reco"][f"{param}[{i}]"] = str(trial.params[key])
        
        with open(OPTIMIZED_CONFIG_FILE, "w") as f:
            config.write(f, space_around_delimiters=False)
        
        logger.info(f"Optimized configuration saved to {OPTIMIZED_CONFIG_FILE}")
        
    except Exception as e:
        logger.error(f"Error saving configuration: {e}")
        raise


def print_study_summary(study: optuna.Study) -> None:
    """Print a summary of the study"""
    print("\n" + "="*60)
    print("OPTIMIZATION RESULTS")
    print("="*60)
    print(f"\nBest trial: #{study.best_trial.number}")
    print(f"Metric value: {study.best_trial.value:.6f}")
    print("\nBest parameters:")
    for key, value in study.best_trial.params.items():
        print(f"  {key:30s} = {value}")
    print("\n" + "="*60 + "\n")


def main():
    """Main function"""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Optimize MS Track Finder CA parameters using Optuna"
    )
    parser.add_argument(
        "-f",
        "--fresh",
        action="store_true",
        help="Start a fresh optimization study (delete existing database)"
    )

    parser.add_argument(
        "-n",
        "--n_trials",
        type=int,
        default=10,
        help="Number of trials to run (default: 10)"
    )
    parser.add_argument(
        "--nevents",
        type=int,
        default=10000,
        help="Number of events to process in each trial (default: 1000)"
    )
    args = parser.parse_args()

    n_trials = args.n_trials  # Increased for better exploration
    n_events = args.nevents

    logger.info("="*60)
    logger.info("Starting MS Track Finder CA optimization")
    logger.info("="*60)
    
    # Handle fresh start
    if args.fresh:
        db_file = Path("optuna_ms_ca.db")
        if db_file.exists():
            db_file.unlink()
            logger.info("Deleted existing database for fresh start")
    
    # Validate required files
    if not validate_files():
        logger.error("File validation failed. Exiting.")
        return
    
    try:
        # Read first trial from INI file
        first_trial = get_first_trial_from_ini()
        logger.info("First trial parameters:")
        for key, value in first_trial.items():
            logger.info(f"  {key}: {value}")
        
        # Create or load Optuna study with persistence
        study = optuna.create_study(
            storage=DB_PATH,
            study_name=STUDY_NAME,
            load_if_exists=True,
            direction="maximize"
        )
        
        # Enqueue first trial with current parameters
        if len(study.trials) == 0:
            study.enqueue_trial(first_trial)
            logger.info("First trial enqueued")
        else:
            logger.info(f"Resuming existing study with {len(study.trials)} completed trials")
        
        # Run optimization
        logger.info(f"Starting optimization with {n_trials} trials...")
        
        study.optimize(
            lambda trial: objective(trial, n_events), 
            n_trials=n_trials,
            show_progress_bar=True,
            catch=(Exception,)  # Continue even if some trials fail
        )
        
        # Print results
        print_study_summary(study)
        
        # Save best configuration
        save_best_config(study)
        
        logger.info("Optimization completed successfully!")
        
    except KeyboardInterrupt:
        logger.warning("Optimization interrupted by user")
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        raise


if __name__ == "__main__":
    main()