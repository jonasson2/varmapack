#ifndef TESTS_H
#define TESTS_H

#include <stdbool.h>

#define TEST_NREP_DEFAULT 10000
extern int TEST_NREP;

void TestTestcase(void);
void Testvarmapack_specrad(void);
bool TestAgainstMatlab(void);
void TestFindCFindG(void);
void TestPsi(void);
void TestAcvf(void);
void TestAutocov(void);
void TestAutocovEdgeCases(void);
void TestCov2corr(void);
void TestSimEdgeCases(void);
void TestSimxEdgeCases(void);
void TestFromMatlab(void);
void TestPsdCondCov(void);
void TestTestcasex(void);
void TestLyapunov(void);
void TestDistribution(void);
#endif /* TESTS_H */
