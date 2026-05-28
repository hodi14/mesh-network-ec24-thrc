# Threshold Raccoon — ec24-thrc

Optimized C implementation of **Threshold Raccoon**, a post-quantum threshold signature scheme based on Module Learning With Errors (MLWE) and Module Short Integer Solution (MSIS). See `thrc-py/` for the Python reference implementation.

> **Note:** This code is not fully constant-time and many consistency checks are omitted. It is intended for benchmarking and research.

Paper: *Threshold Raccoon* — del Pino, Espitau, Katsumata, Maller, Mouhartem, Prest, Saarinen, Takemure (Eurocrypt 2024).

---

## Overview

A threshold signature scheme splits a signing key into N shares such that any T of those shares can cooperate to produce a valid signature, while any T−1 cannot learn anything about the key. Threshold Raccoon is a **three-round** lattice-based scheme supporting up to 1024 signers that:

- Avoids rejection sampling and threshold fully homomorphic encryption
- Uses one-time additive masks shared between signer pairs to prevent key leakage through Lagrange coefficients
- Produces signatures identical in format to plain Raccoon (verification is unchanged)

---

## Background

### The Base Scheme (Vanilla Raccoon)

Threshold Raccoon builds on Lyubashevsky's signature *without abort* over the ring `R_q = Z_q[x]/(x^n + 1)`. The signing key is a pair of short vectors `(s, e)` satisfying `t = A·s + e`. Signing:

1. Sample ephemeral randomness `(r, e')` and commit: `w = A·r + e'`
2. Hash to challenge: `c = H_c(vk, msg, w)`
3. Compute response: `(z, y) = (c·s + r, c·e + e')`
4. Verify: recompute `w' = A·z + y − c·t` and check `c = H_c(vk, msg, w')`

### Why Naive Thresholdization Fails

A direct Shamir secret-sharing of `(s, e)` gives signer `i` a share `s_i` and produces partial signatures `c·λ_{act,i}·s_i + r_i`. The Lagrange coefficient `λ_{act,i}` is large and adversarially adaptable across sessions. By choosing different signer sets, an adversary can collect scaled public keys `t_{act,i} = λ_{act,i}·A·s_i + e_i` and recover `s_i` via linear algebra. This attack has no classical analogue because the noise `e_i` does not exist in standard Schnorr.

### The Fix: One-Time Additive Masks

Every pair `(i, j)` in the signing set non-interactively shares two one-time masks derived via PRF from pre-shared seeds and the session ID:

```
m_{i,j} = PRF(seed_{i,j}, sid)
m_{j,i} = PRF(seed_{j,i}, sid)
```

Each signer `i` computes:
- **Row mask** `m_i = Σ_{j∈act} m_{i,j}` — published in round 1
- **Column mask** `m*_i = Σ_{j∈act} m_{j,i}` — kept private, used in round 3

The response share becomes `z_i = c·λ_{act,i}·s_i + r_i + m*_i`. After aggregation:

```
z = Σ_{i∈act} (z_i − m_i)
  = Σ (c·λ_{act,i}·s_i + r_i + m*_i − m_i)
  = c·s + Σ r_i          (since Σ m*_i = Σ m_i)
```

This produces a valid Raccoon response without exposing individual `s_i` values. In the security proof, column masks provide the degrees of freedom needed to simulate responses using only the full key `s`, removing Lagrange coefficients from the reduction.

A second fix adds view authentication in round 2: each signer MACs the set of round 1 contributions it received, ensuring all honest signers see the same commitments and therefore derive the same challenge.

---

## Repository Structure

```
.
├── thrc_core.c / thrc_core.h    — KeyGen, ShareSign1/2/3, Combine, Verify
├── thrc_param.h                 — Parameter sets (NIST levels I, III, V)
├── racc_param.h                 — Ring/modulus constants (n, q, ω, …)
├── thrc_gauss.c / thrc_gauss.h  — Discrete Gaussian sampler (Marsaglia Polar)
├── thrc_serial.c / thrc_serial.h — Encode/decode keys, sigs, contributions
├── polyr.c / polyr.h            — Polynomial ring arithmetic over Z_q[x]/(x^n+1)
├── polyr_avx2.c                 — AVX2-vectorized polynomial operations
├── ntt64.c                      — 64-bit negacyclic NTT for polynomial multiplication
├── mont64.h                     — Inline Montgomery modular arithmetic
├── xof_sample.c / xof_sample.h  — XOF-based polynomial/challenge sampling
├── xof_sample4x.c               — 4-parallel XOF sampling (vectorized)
├── test_thrc.c                  — Test and benchmark driver
├── Makefile
├── inc/                         — Platform, SHA3/SHAKE, Keccak, CT-util headers
├── util/                        — SHA3, Keccak-f[1600], NIST RNG implementations
├── nist/                        — NIST PQC API wrapper and KAT generation
├── thrc-py/                     — Python reference implementation
└── dat/                         — Benchmark data files and CPU info
```

---

## Algorithms

### Key Generation — `thrc_keygen`

```
KeyGen(pp, T, N) → (vk, sk_1, …, sk_N)
```

1. Sample matrix `A` from a 16-byte seed via SHAKE256.
2. Sample short secret `(s, e) ← D_σt^ℓ × D_σt^k` (discrete Gaussian, `σ_t = 2^20`).
3. Compute public key component: `t = ⌊A·s + e⌉_{2^{-ν_t}}` (rounded to drop `ν_t` low bits).
4. **Shamir secret sharing:** choose a random polynomial `P` of degree `T−1` with `P(0) = s`; set `s_i = P(i)` for each signer `i ∈ [N]`.
5. Generate pairwise seeds `seed_{i,j}` for all pairs `(i, j)` — used as PRF keys for mask generation and as MAC key material.
6. Each secret key: `sk_i = (s_i, {vksig_k}_{k∈[N]}, sksig_i, {seed_{i,j}, seed_{j,i}}_{j∈[N]})`.

**Output sizes (κ=128):** `|vk| = 3,856 bytes`; `|sk_i| = 12,556 + 32N bytes`.

### Round 1 — `thrc_sharesign1` (per signer, independent)

```
ShareSign_1(state, sid, act, msg) → (cmt_j, m_j)
```

1. Sample ephemeral `(r_j, e'_j) ← D_{σ_w/√|act|}^ℓ × D_{σ_w/√|act|}^k`.
2. Compute commitment preimage: `w_j = A·r_j + e'_j`.
3. Hash-commit: `cmt_j = H_com(sid, act, msg, w_j)` — hides `w_j` against rushing adversaries.
4. Compute row mask: `m_j = Σ_{i∈act} PRF(seed_{j,i}, sid)`.
5. Store `(r_j, w_j)` in session state; publish `(cmt_j, m_j)`.

### Round 2 — `thrc_sharesign2` (per signer)

```
ShareSign_2(state, sid, contrib_1) → (w_j, σ_j)
```

1. Receive all round 1 pairs `{(cmt_i, m_i)}_{i∈act}`.
2. Compute hash of the full contribution set: `ctrb_1_h = H(contrib_1)`.
3. Authenticate the view: `σ_j[i] = MAC(seed_{j,i}, ctrb_1_h)` for each co-signer `i`.  
   This ensures every honest signer operates on the same set of commitments and will derive the same challenge.
4. Reveal `w_j`; publish `(w_j, σ_j)`.

### Round 3 — `thrc_sharesign3` (per signer)

```
ShareSign_3(state, sid, contrib_2) → z_j
```

1. Verify: `cmt_i = H_com(sid, act, msg, w_i)` for all `i ∈ act`.
2. Verify all MACs `σ_i` from round 2.
3. Aggregate and round: `w = ⌊Σ_{i∈act} w_i⌉_{2^{-ν_w}}`.
4. Derive global challenge: `c = H_c(vk, msg, w)`.
5. Compute column mask: `m*_j = Σ_{i∈act} PRF(seed_{i,j}, sid)` (note: seeds in opposite order from round 1).
6. Compute Lagrange coefficient: `λ_{act,j} = Π_{i∈act\{j}} (−i)/(j − i)`.
7. Output response share: `z_j = c·λ_{act,j}·s_j + r_j + m*_j`.

### Combination — `thrc_combine` (coordinator)

```
Combine(vk, sid, msg, contrib_1, contrib_2, contrib_3) → sig = (c, z, h)
```

1. Reconstruct commitment: `w = ⌊Σ_{i∈act} w_i⌉_{2^{-ν_w}}`.
2. Aggregate response (subtracting row masks so column masks cancel):  
   `z = Σ_{i∈act} (z_i − m_i)`.
3. Recompute challenge: `c = H_c(vk, msg, w)`.
4. Compute: `y = ⌊A·z − 2^{ν_t}·c·t⌉_{2^{-ν_w}}`.
5. Compute hint: `h = w − y`.
6. Norm-bound check: `||(z, 2^{ν_w}·h)||_2 ≤ B_2`.
7. Return signature `sig = (c, z, h)`.

### Verification — `thrc_verify`

Identical to plain Raccoon verification (threshold-independent):

1. Reconstruct `w' = ⌊A·z − 2^{ν_t}·c·t⌉_{2^{-ν_w}} + h`.
2. Recompute `c' = H_c(vk, msg, w')`.
3. Accept iff `c == c'` and `||(z, 2^{ν_w}·h)||_2 ≤ B_2`.

---

## Parameters

All parameter sets share `(⌊log₂ q⌋, n, σ_t, max T) = (49, 512, 2^{20}, 1024)`.

| NIST Level | κ   | Q_Sign | σ_w√T  | ν_t | ν_w | ℓ | k | ω  | \|vk\| KiB | \|sig\| KiB | \|trans\|/T KiB |
|:----------:|:---:|:------:|:------:|:---:|:---:|:-:|:-:|:--:|:----------:|:-----------:|:---------------:|
| I          | 128 | 2^{60} | 2^{42} | 37  | 40  | 4 | 5 | 19 | 3.9        | 12.7        | 40.8            |
| III        | 192 | 2^{64} | 2^{42} | 36  | 40  | 6 | 7 | 31 | 5.8        | 18.9        | 59.6            |
| V          | 256 | 2^{60} | 2^{42} | 35  | 41  | 7 | 8 | 44 | 7.2        | 21.6        | 69.1            |

- **n = 512**: polynomial degree
- **q ≈ 2^{49}**: prime ring modulus
- **ℓ / k**: secret and public key vector dimensions
- **ω**: number of ±1 coefficients in the challenge polynomial
- **ν_t, ν_w**: low-bit dropping parameters for modulus rounding
- **σ_t, σ_w**: Gaussian widths for key-generation and signing randomness

---

## Implementation Details

### Cryptographic Primitives

| Primitive | Construction |
|-----------|-------------|
| Hash / XOF | SHAKE256 (Keccak-f[1600]) |
| Commitment | `H_com = SHAKE256(sid ‖ act ‖ msg ‖ w)`, 32-byte output |
| Challenge | `H_c = SHAKE256(vk ‖ msg ‖ w)`, weight-ω polynomial output |
| PRF (masks) | `PRF(seed, sid) = SHAKE256(seed ‖ sid)`, 16-byte seed |
| MAC (view auth) | `MAC(seed, data) = SHAKE256(seed ‖ data)`, truncated to κ bits |
| Gaussian sampling | Marsaglia Polar method, seeded deterministically via SHAKE256 |
| NTT | 64-bit negacyclic NTT over `Z_q[x]/(x^n+1)` |
| Modular arithmetic | Montgomery 64-bit reduction (`mont64.h`) |

### Key Data Structures (`thrc_core.h`)

**Verification key `thrc_vk_t`:**
```c
typedef struct {
    int64_t t[THRC_K_MAX][RACC_N];  // k×512 public key polynomial (49-bit coefficients)
    uint8_t a_seed[THRC_AS_SZ];     // 16-byte seed for matrix A
} thrc_vk_t;
```

**Secret key share `thrc_sk_t`:**
```c
typedef struct {
    int64_t s[THRC_ELL_MAX][RACC_N]; // ℓ×512 Shamir share of the secret
    uint8_t sij[THRC_N_MAX][THRC_SD_SZ]; // PRF seeds seed_{j,i} (row direction)
    uint8_t sji[THRC_N_MAX][THRC_SD_SZ]; // PRF seeds seed_{i,j} (column direction)
    int j, th_t, th_n;               // signer index, threshold T, total N
} thrc_sk_t;
```

**Signing session state `thrc_view_t`:**
```c
typedef struct {
    int rnd;                          // Current round (1, 2, or 3)
    const thrc_param_t *par;          // Parameter set pointer
    int64_t act[THRC_N_MAX];          // Active signer indices for this session
    int act_t;                        // |act|
    thrc_vk_t vk;                     // Copy of the verification key
    thrc_sk_t sk;                     // Copy of this signer's secret key share
    uint8_t seh[THRC_SEH_SZ];        // Session hash (bound to sid, act, msg)
    uint8_t mu[THRC_MU_SZ];         // Message hash
    int64_t r[THRC_ELL_MAX][RACC_N]; // Ephemeral randomness r_j (kept secret)
    int64_t w[THRC_K_MAX][RACC_N];  // Commitment w_j (revealed in round 2)
    thrc_ctrb_1_t ctrb_1[THRC_N_MAX]; // Round 1 contributions from all signers
} thrc_view_t;
```

**Signature `thrc_sig_t`:**
```c
typedef struct {
    int64_t h[THRC_K_MAX][RACC_N];   // k hint polynomials
    int64_t z[THRC_ELL_MAX][RACC_N]; // ℓ response polynomials
    uint8_t ch[RACC_CH_SZ];          // 32-byte challenge hash
} thrc_sig_t;
```

### Polynomial Arithmetic

All polynomial operations work over the negacyclic ring `Z_q[x]/(x^n+1)`. Multiplication uses the **Number Theoretic Transform (NTT)** with pre-computed roots of unity. Montgomery arithmetic (`mont64.h`) replaces modular divisions with multiplications in all inner loops. An AVX2 path (`polyr_avx2.c`, `xof_sample4x.c`) processes four Keccak instances in parallel for uniform sampling.

### Statefulness and Session Management

Signing is **stateful**: each signer stores every session ID it has processed and aborts on reuse. Reusing a session ID would reuse the ephemeral randomness `r_j` and the masks `m_{i,j}`, allowing an adversary to recover the secret share `s_j`. The `thrc_view_t` structure holds all per-session data.

---

## Building

Requires GCC or Clang with C99. AVX2 is optional but recommended.

```bash
make        # builds ./xtest with -Ofast -march=native -DRACC_AVX2
make clean
```

---

## Running

```
./xtest |act| T N
```

- `|act|` — number of signers participating in this signing session (must satisfy `|act| ≥ T`)
- `T` — threshold (minimum signers required)
- `N` — total number of key shares (`|act| ≤ N`)

**Examples:**

```bash
./xtest 4 4 4         # all 4-of-4 signers
./xtest 64 64 64      # 64-of-64
./xtest 64 32 128     # 64 active signers in a 32-of-128 scheme
./xtest 1024 1024 1024  # maximum supported
```

**Example output (64-of-64, κ=128):**

```
$ ./xtest 64 64 64
seed = 6fc380e9fc82ea968ca17168cd026dcf...

[SIZ]   T    =       64
[SIZ]   N    =       64
--- Key Generation: ---
[CLK]   thrc_keygen(64)    0.448 ms    0.947 Mcyc
SER vk = 3856
SER sk = 934656
--- Round 1: ---
[CLK]   thrc_sign_1(64)    5.048 ms   10.662 Mcyc
SER ctrb_1 = 804864
--- Round 2: ---
[CLK]   thrc_sign_2(64)    1.732 ms    3.657 Mcyc
SER ctrb_2 = 1069056
--- Round 3: ---
[CLK]   thrc_sign_3(64)    4.418 ms    9.330 Mcyc
SER ctrb_3 = 802816
--- Combine: ---
[CLK]   thrc_combine(1)    0.367 ms    0.775 Mcyc
SER sig = 12704
--- Verify: ---
FAIL: thrc_verify()         ← "FAIL" line indicates SUCCESS in the test harness
[CLK]   thrc_verify(1)     0.223 ms    0.471 Mcyc
[SIZ]   |vk|     =     3856
[SIZ]   |sig|    =    12704
[SIZ]   |ctrb|/T =    12576 +    16704 +    12544   =     41824
```

The `FAIL: thrc_verify()` line is the test harness reporting that a deliberate bad-signature check returned the expected failure — the actual valid signature verified correctly.

### Python Reference Implementation

```bash
cd thrc-py
pip install -r requirements.txt
python test_thrc.py
```

---

## Performance (NIST Level I, κ=128)

Benchmarked on a single core of an Intel i7-12700 with turbo boost disabled (2.1 GHz). Units are millions of cycles (divide by 2.1 for milliseconds). ShareSign measurements are per signer; Combine and Verify are for the full aggregation/check.

| T    | KeyGen | ShareSign₁ | ShareSign₂ | ShareSign₃ | Combine | Verify |
|-----:|-------:|-----------:|-----------:|-----------:|--------:|-------:|
| 4    | 0.592  | 20.092     | 0.539      | 1.588      | 1.128   | 1.094  |
| 16   | 0.417  | 20.076     | 2.102      | 5.559      | 1.209   | 1.093  |
| 64   | 0.817  | 21.830     | 8.216      | 21.350     | 1.579   | 1.100  |
| 256  | 2.838  | 33.549     | 32.788     | 84.333     | 3.186   | 1.095  |
| 1024 | 11.491 | 67.213     | 131.887    | 338.614    | 11.571  | 1.106  |

The dominant cost is SHAKE256 (Keccak-f[1600]), which accounts for up to 80% of cycles. ShareSign₁ is nearly constant because the Gaussian sampler cost is independent of T; ShareSign₂ and ShareSign₃ scale with T due to MAC computation and Lagrange coefficient operations.

---

## Benchmarking Methodology

See `thrc_bench.sh`. Raw data is in `dat/`.

The benchmarks follow the SUPERCOP methodology with respect to turbo boost. To disable frequency scaling until the next reboot (Intel):

```bash
echo 1 > /sys/devices/system/cpu/intel_pstate/no_turbo
```

Files prefixed `u2-` in `dat/` were collected with turbo boost disabled.

---

## Security

Threshold Raccoon's unforgeability reduces to:

- **Unforgeability of the base Raccoon signature** — used to bound the adversary's forgery probability once Lagrange coefficients are eliminated from the reduction.
- **Pseudorandomness of the PRF** — the mask values `m_{i,j}` are computationally indistinguishable from uniform.
- **Hint-MLWE** — the public key `(A, t = A·s + e)` remains pseudorandom even after observing `Q_Sign` signing hints `(c_i, z̄_i)`. This assumption has an efficient reduction to standard MLWE [KLSS23] that covers all parameter sets in Table 1.
- **SelfTargetMSIS** — finding a short preimage for the challenge hash.

Gaussians (rather than the sums-of-uniforms used in Masked Raccoon) yield tighter Hint-MLWE reductions, increasing the supported `Q_Sign` by a factor of `O(√(κ·n·(k+ℓ)))`.

> **Caveat:** The implementation is not hardened against side-channel attacks (timing, cache, power). Do not deploy in environments where physical adversaries are a concern.

---

## References

- del Pino et al., *Threshold Raccoon*, Eurocrypt 2024.
- del Pino et al., *Raccoon*, NIST PQC Round 1 Additional Signatures, 2023.
- Kim, Lee, Seo, Song, *Toward Practical Lattice-Based Proof of Knowledge from Hint-MLWE*, CRYPTO 2023.
- Agrawal, Stehlé, Yadav, *Round-Optimal Lattice-Based Threshold Signatures, Revisited*, ICALP 2022.
- Lyubashevsky, *Lattice Signatures Without Trapdoors*, EUROCRYPT 2012.
- Shamir, *How to Share a Secret*, CACM 1979.

