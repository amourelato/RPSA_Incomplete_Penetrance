from pathlib import Path

def get_workflow_dir():
    """
    Returns the absolute Path to the workflow1 directory,
    based on the location of this file.

    Usage:
        In any script inside RPSA_eQTL_Fine_Mapping/scripts/, call this to get base path.
    """
    script_dir = Path(__file__).resolve().parent
    workflow_dir = script_dir.parent 
    return workflow_dir