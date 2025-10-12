test_that("OmegaBuild matches block constructions for r=2, p=3, n=7, q in {2,3,4}", {
  r <- 2; p <- 3; n <- 7

  S0 <- matrix(c(3,8,4,8), 2, 2)  # be careful: R fills by column; values taken row-wise in your MATLAB
  S1 <- matrix(c(3,1,4,1), 2, 2)
  S2 <- matrix(c(3,2,4,2), 2, 2)
  S3 <- matrix(c(3,3,4,3), 2, 2)

  G0 <- matrix(c(4,8,5,8), 2, 2)
  G1 <- matrix(c(4,1,5,1), 2, 2)
  G2 <- matrix(c(4,2,5,2), 2, 2)
  G3 <- matrix(c(4,3,5,3), 2, 2)
  G4 <- matrix(c(4,4,5,4), 2, 2)

  W0 <- matrix(c(5,8,6,8), 2, 2)
  W1 <- matrix(c(5,1,6,1), 2, 2)
  W2 <- matrix(c(5,2,6,2), 2, 2)
  W3 <- matrix(c(5,3,6,3), 2, 2)
  W4 <- matrix(c(5,4,6,4), 2, 2)

  O <- matrix(0, r, r)

  # Convenience constructors for S, G, W concatenation (r x r*(k+1))
  Scat <- cbind(S0, S1, S2, S3)
  Gcat_2 <- cbind(G0, G1, G2)
  Wcat_2 <- cbind(W0, W1, W2)
  Gcat_4 <- cbind(G0, G1, G2, G3, G4)
  Wcat_4 <- cbind(W0, W1, W2, W3, W4)

  ## Case 1: q < p (q=2)
  q <- 2
  out <- omega_build(Scat, Gcat_2, Wcat_2, p = p, q = q, r = r, n = n)
  Su  <- out$Su
  Olow <- out$Olow

  SuOK <- rowblk(
            blk(S0, O,  O),
            blk(S1, S0, O),
            blk(S2, S1, S0)
          )
  OlowOK <- rowblk(
              blk(G2, G1, W0),
              blk(G2, W1, W0),
              blk(W2, W1, W0),
              blk(W2, W1, W0)
            )

  expect_equal(tril(Su), tril(SuOK))
  expect_equal(Olow, OlowOK)

  ## Case 2: q = p + 1 (q=4 with p=3)
  q <- 4
  outB <- omega_build(Scat, Gcat_4, Wcat_4, p = p, q = q, r = r, n = n)
  SuB  <- outB$Su
  OlowB <- outB$Olow

  SOKB <- rowblk(
            blk(S0, O,  O,  O),
            blk(S1, S0, O,  O),
            blk(S2, S1, S0, O),
            blk(G3, G2, G1, W0)
          )
  OlowOKB <- rowblk(
               blk(G4, G3, G2, W1, W0),
               blk(G4, G3, W2, W1, W0),
               blk(G4, W3, W2, W1, W0)
             )

  expect_equal(tril(SuB), tril(SOKB))
  expect_equal(OlowB, OlowOKB)

  ## Case 3: q = p + 2 but using SC = {S0,S1,S2} and p=2
  q <- 4
  pC <- 2
  SCat <- cbind(S0, S1, S2)
  outC <- omega_build(SCat, Gcat_4, Wcat_4, p = pC, q = q, r = r, n = n)
  SuC  <- outC$Su
  OlowC <- outC$Olow

  SOKC <- rowblk(
            blk(S0, O,  O,  O),
            blk(S1, S0, O,  O),
            blk(G2, G1, W0, O),
            blk(G3, G2, W1, W0)
          )
  OlowOKC <- rowblk(
               blk(G4, G3, W2, W1, W0),
               blk(G4, W3, W2, W1, W0),
               blk(W4, W3, W2, W1, W0)
             )

  expect_equal(tril(SuC), tril(SOKC))
  expect_equal(OlowC, OlowOKC)
})
