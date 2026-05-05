#ifndef TESTS_H
#define TESTS_H

#include <stdbool.h>

void TestExtraUtil(void);
void TestError(void);
void TestTestcase(void);
void Testvarmapack_specrad(void);
bool TestAgainstMatlab(void);
void TestFindCG(void);
void TestPsi(void);
void TestAcvf(void);
void TestAutocov(void);
void TestAutocovEdgeCases(void);
void TestFromMatlab(void);
void TestPsdCondCov(void);
#endif /* TESTS_H */
