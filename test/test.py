from yaml_utils import *
from crab_manager import crab_job_manager

manager = crab_job_manager()

manager.setup_crab('../test/configs/example_config.yml')
