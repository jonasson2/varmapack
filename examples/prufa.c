#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

double ScalarFromMatlab(char *filename, char *var);

static bool print_scalar(char *filename, char *name)
{
  double value = ScalarFromMatlab(filename, name);
  if (isnan(value)) {
    printf("Could not read %s from %s\n", name, filename);
    return false;
  }
  printf("%s = %.6f\n", name, value);
  return true;
}

int main(int argc, char *argv[])
{
  const char default_file[] = "examples/prufa_data.txt";
  char filename[256];
  bool ok = true;

  snprintf(filename, sizeof filename, "%s", default_file);
  if (argc > 1) {
    snprintf(filename, sizeof filename, "%s", argv[1]);
  }

  printf("Reading scalars from %s\n", filename);

  ok = print_scalar(filename, "alpha") && ok;
  ok = print_scalar(filename, "gamma") && ok;

  return ok ? 0 : 1;
}
