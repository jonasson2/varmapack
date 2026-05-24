#ifndef TESTS_H
#define TESTS_H

#include <stdbool.h>

void TestExtraUtil(void);
void TestTestcase(void);
void Testvarmapack_specrad(void);
bool TestAgainstMatlab(void);
void TestFindCFindG(void);
void TestPsi(void);
void TestAcvf(void);
void TestAutocov(void);
void TestAutocovEdgeCases(void);
void TestSimEdgeCases(void);
void TestSimxEdgeCases(void);
void TestFromMatlab(void);
void TestPsdCondCov(void);
#endif /* TESTS_H */
