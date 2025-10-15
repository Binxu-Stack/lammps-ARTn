# Install/unInstall package classes in LAMMPS

if (test $1 = 1) then

  cp -p min_artn.h ..
  cp -p min_artn.cpp ..

elif (test $1 = 0) then

  rm ../min_artn.h
  rm ../min_artn.cpp

fi
