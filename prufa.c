int main(void) {
  int A[] = {0,0,0};
  A[1] = (static int[4]) {0, 1, 2, 3};
  return A[1][0];
}
