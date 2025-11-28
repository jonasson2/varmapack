#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include "../tests/FromMatlab.h"

static bool print_scalar(FILE *f, const char *filename, char *name) {
  double value = 0.0;
  bool ok = DoubleFromMatlab(f, name, &value);
  if (!ok) {
    printf("Could not read %s from %s\n", name, filename);
    return false;
  }
  printf("%s = %.6f\n", name, value);
  return true;
}

int main(int argc, char *argv[]) {
  const char default_file[] = "prufa_data.txt";
  char filename[256];
  bool ok = true;

  snprintf(filename, sizeof filename, "%s", default_file);
  if (argc > 1) {
    snprintf(filename, sizeof filename, "%s", argv[1]);
  }

  printf("Reading scalars from %s\n", filename);

  FILE *f = fopen(filename, "r");
  if (!f) {
    printf("Could not open %s\n", filename);
    return 1;
  }

  ok = print_scalar(f, filename, "alpha") && ok;
  ok = print_scalar(f, filename, "gamma") && ok;

  fclose(f);

  return ok ? 0 : 1;
}
