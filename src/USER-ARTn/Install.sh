# Install/unInstall package classes in LAMMPS

# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

action () {
  if (test $mode = 0) then
    rm -f ../$1
  fi
}

# all package files with no dependencies

for file in *.cpp *.h; do
  test -f ${file} && action $file
done

