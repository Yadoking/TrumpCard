#!/bin/bash

if [ ! -d Delphes ]; then
  #ln -s ../Delphes-3.4.1 Delphes
  ln -s /home/jhgoh/work/Delphes/Delphes-3.4.1 Delphes
fi

echo $LD_LIBRARY_PATH | grep -q Delphes
if [ $? -ne 0 ]; then
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`readlink -f Delphes`
  export PATH=$PATH:`readlink -f Delphes`
fi


