#! /bin/sh

# Basically, we're going to do into each directory and run all
# codes that are of the pattern 'test*.py'

# Copy the test directory to /tmp and run the tests
home="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

tmp_dir="/tmp/unittests_$$"
mkdir $tmp_dir
cp -r $home/*_tests $tmp_dir
cd $tmp_dir

# In case a test needs to get back to the source directory
export STEAM_TEST_HOME=$home

pass=0
fail=0

dirs=`ls -d *_tests`
for dir in $dirs
do
   echo " "
   echo " "
   echo " "
   echo "============================================================"
   echo "   $dir"
   echo "============================================================"
   echo " "
   cd $dir
   # https://docs.python.org/3/library/unittest.html#command-line-interface
   python -m unittest discover -vb | tee output.txt
   if [ ${PIPESTATUS[0]} -eq 0 ]
   then
   #  echo "recording a pass... ${PIPESTATUS[0]}"
      ((pass++))
   else
   #  echo "recording a fail... ${PIPESTATUS[0]}"
      ((fail++))
   fi   
   cd -
   echo " "
   echo " "
done

echo " "
echo "SUMMARY:"
echo " - Pass: $pass"
echo " - Fail: $fail"
echo " "
echo " Details are in $tmp_dir/*/output.txt"
echo " "

if [ $fail != 0 ] ; then
   exit $fail
fi
