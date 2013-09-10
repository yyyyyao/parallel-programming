#for N in 1024 2048 4096 8192; do
for N in 64 128 256 512; do
  echo "---------------"
  echo N:$N
  COMPUTE_PROFILE='1' COMPUTE_PROFILE_CSV='1' COMPUTE_PROFILE_CONFIG='event.prof' ./main $N
  sed -n 6p ./opencl_profile_0.log
  sed -n 9p ./opencl_profile_0.log
done
