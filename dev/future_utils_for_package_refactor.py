# DEPENDENCY CHECK WILL PROBABLY BE USEFUL IF THIS WAS NOT STANDALONE
import shutil
# from apps.binning.unitem.unitem import UniteM
# from apps.read_mapping.coverm.coverm import CoverM

PATH_PKGS = ['metabat', 'run_MaxBin.pl', "UniteM",
             "coverm"]


def check_dependencies(ref_pkgs: List[str]) -> None:
    """
    Check if all dependencies are on system path and are executable
    """
    acc = [shutil.which(pkg) for pkg in ref_pkgs if shutil.which(pkg) is not None]
    for missing_pkg in acc:
        print(f'{missing_pkg} is not installed')
    return None

if __name__ == "__main__":
    pass