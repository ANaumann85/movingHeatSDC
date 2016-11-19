for d in $(ls -d kiter_*); do
  cd $d
  python ../accuracy.py
  cd ..
done
