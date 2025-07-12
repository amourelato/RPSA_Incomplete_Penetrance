from pathlib import Path

def get_workflow_dir():
    """
    Returns the absolute Path to the WT_RPSA_Analysis directory,
    based on the location of this file.

    Usage:
        In any script inside WT_RPSA_Analysis/scripts/, call this to get base path.
    """
    
    script_dir = Path(__file__).resolve().parent
    workflow_dir = script_dir.parent  
    return workflow_dir