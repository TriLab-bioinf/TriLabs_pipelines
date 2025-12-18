def get_chr_ids(intervals: str) -> list:
    """
    Generate a list of chromosome ids from the intervals file.
    """
    chrIDs: list = []
    with open(intervals, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                chrIDs.append(line.split()[0])
    chrIDs = list(set(chrIDs))  # Remove duplicates
    return chrIDs


def process_db_action(_db_action: str) -> str:
    action: str = ""
    if _db_action == "create":
        action = "--genomicsdb-workspace-path"
    elif _db_action == "update":
        action = "--genomicsdb-update-workspace-path"
    else:
        raise ValueError(
            "invalid option provided to 'params.db_action'; please choose either 'create' or 'update'."
        )
    
    return action

def fetch_genomeDB_id(db_version_file: str) -> str:
    db_id: str = '' 
    with open(db_version_file, 'r') as f:
        db_id = f.readline.rstrip()
    
    return db_id