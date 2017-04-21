echo "Copying files.."
cp $cuode/Makefile .
cp $cuode/input.h .
echo ".. copied Makefile and input.h"
echo "Now linking .." 
mkdir -p src
cp $cuode/Makefile .
ln -sf $cuode/src/Makefile  src/
ln -sf $cuode/src/solve.cu  src/
ln -sf $cuode/src/ode.h   src/
ln -sf $cuode/src/model.h   src/
ln -sf $cuode/src/ode.cu    src/
ln -sf $cuode/src/CUDA    src/ 
ln -sf $cuode/src/models    src/ 
ln -sf $cuode/src/modules    src/ 
ln -sf $cuode/src/modules    src/ 
ln -sf $cuode/src/toolbox    src/ 
echo "...Linking done"
echo "Choose models and parameters and run make"
