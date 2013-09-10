__kernel void vector_add(__global const int *A,
                                  __global const int *B,
                                  __global int *C)
{
  int i = get_global_id(0);
  int sum = 0;

  if(A[i] % 2) {
    sum = 5;
    //if(A[i] / 7 == 0) sum += 3;
  } else {
    sum = 0;
  }
  C[i] = A[i] + B[i] + sum;
}
