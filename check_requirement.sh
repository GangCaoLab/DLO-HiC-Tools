#!/bin/bash

unsatisfy=0

# check python packages
python -c "import regex"      &> /dev/null  || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] Python package \"regex\" uninstall"; }
python -c "import HTSeq"      &> /dev/null  || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] Python package \"HTSeq\" uninstall, see \"http://www-huber.embl.de/HTSeq/doc/overview.html\""; }
python -c "import numpy"      &> /dev/null  || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] Python package \"regex\" uninstall"; }
python -c "import matplotlib" &> /dev/null  || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] Python package \"matplotlib\" uninstall"; }
python -c "import iced"       &> /dev/null  || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] Python package \"hiclib/iced\" uninstall, see \"https://github.com/hiclib/iced/blob/master/doc/tutorial/index.rst\"."; }

# check other programs
which fastx_clipper &> /dev/null || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] fastx_tools uninstall, see \"http://hannonlab.cshl.edu/fastx_toolkit/download.html\""; }
which bwa           &> /dev/null || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] bwa uninstall"; }
which bedtools      &> /dev/null || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] bedtools uninstall"; }
which samtools      &> /dev/null || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] samtools uninstall"; }
which parallel      &> /dev/null || { unsatisfy=$(($unsatisfy+1)) ; echo "[warning] parallel uninstall"; }

if [ $unsatisfy -gt 0 ]; then
    echo "[warning] $unsatisfy softwares/packages uninstall!"
fi

exit $unsatisfy
