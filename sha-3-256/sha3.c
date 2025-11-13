/*
 * SHA3-256 file hasher (Keccak-f[1600], NIST SHA-3, padding 0x06..0x80).
 *
 * Key points:
 * - Sponge parameters: rate=1088 bits (136 bytes), capacity=512 bits (64 bytes), state=1600 bits (200 bytes).
 * - Byte order: lanes (64-bit) are assembled/disassembled in little-endian regardless of host endianness.
 * - Padding: domain separator 0x06 and final bit 0x80 per FIPS 202 (SHA-3), not Keccak-256.
 * - I/O: reads files in binary mode ("rb") in RATE-sized chunks; errors are returned via int status codes.
 * - Security: uses fixed-size local buffers only (state=200B, buffer=136B).
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>

/*
 * Constants defining the Keccak state and SHA3-256 parameters.
 * STATE_SIZE: 1600 bits = 200 bytes total state (5x5 lanes of 64 bits).
 * RATE: bytes absorbed/squeezed per permutation for SHA3-256 (1088/8 = 136).
 * CAPACITY: remaining bytes (512/8 = 64); not used directly (STATE_SIZE = RATE + CAPACITY).
 * DIGEST_SIZE: output size (256 bits = 32 bytes).
 */
#define STATE_SIZE 200
#define RATE 136
#define DIGEST_SIZE 32

/* Keccak-f[1600] round constants for iota step (24 rounds). */
const uint64_t round_constants[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL,
    0x800000000000808aULL, 0x8000000080008000ULL,
    0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008aULL, 0x0000000000000088ULL,
    0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL,
    0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL,
    0x8000000080008081ULL, 0x8000000000008080ULL,
    0x0000000080000001ULL, 0x8000000080008008ULL
};

/* Keccak state as 5x5 array of 64-bit lanes. */
typedef uint64_t state_t[5][5];

/*
 * Left-rotate 64-bit word by n (0 <= n < 64).
 * Uses unsigned shifts; right shift inserts zeros. No UB since n <= 62 in rho.
 */
static inline uint64_t rol(uint64_t x, int n) {
    return (x << n) | (x >> (64 - n));
}

/* Theta: mix columns via parity and rotated neighbors. */
void theta(state_t state) {
    uint64_t C[5], D[5];
    for (int x = 0; x < 5; x++) {
        C[x] = state[x][0] ^ state[x][1] ^ state[x][2] ^ state[x][3] ^ state[x][4];
    }
    for (int x = 0; x < 5; x++) {
        D[x] = C[(x + 4) % 5] ^ rol(C[(x + 1) % 5], 1);
    }
    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            state[x][y] ^= D[x];
        }
    }
}

/* Rho: rotate each lane by a fixed offset. */
void rho(state_t state) {
    static const int rho_offsets[5][5] = {
        {0, 36, 3, 41, 18},
        {1, 44, 10, 45, 2},
        {62, 6, 43, 15, 61},
        {28, 55, 25, 21, 56},
        {27, 20, 39, 8, 14}
    };
    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            state[x][y] = rol(state[x][y], rho_offsets[x][y]);
        }
    }
}

/* Pi: permute lane positions. */
void pi(state_t state) {
    state_t temp;
    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            temp[x][y] = state[(x + 3 * y) % 5][x];
        }
    }
    memcpy(state, temp, sizeof(temp));
}

/* Chi: non-linear step across rows. */
void chi(state_t state) {
    state_t temp;
    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            temp[x][y] = state[x][y] ^ ((~state[(x + 1) % 5][y]) & state[(x + 2) % 5][y]);
        }
    }
    memcpy(state, temp, sizeof(temp));
}

/* Iota: XOR round constant into state[0][0]. */
void iota(state_t state, int round) {
    state[0][0] ^= round_constants[round];
}

/*
 * keccak_f: apply Keccak-f[1600] permutation to the 200-byte state buffer.
 *
 * - Input/output state is a 200-byte array (5x5 lanes x 8 bytes), lanes encoded LE.
 * - Converts bytes -> lanes (LE), performs 24 rounds, writes lanes -> bytes (LE).
 */
void keccak_f(uint8_t *state_bytes) {
    state_t state;
    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            state[x][y] = 0;
            for (int z = 0; z < 8; z++) {
                state[x][y] |= ((uint64_t)state_bytes[(y * 5 + x) * 8 + z]) << (z * 8);
            }
        }
    }

    for (int round = 0; round < 24; round++) {
        theta(state);
        rho(state);
        pi(state);
        chi(state);
        iota(state, round);
    }

    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            for (int z = 0; z < 8; z++) {
                state_bytes[(y * 5 + x) * 8 + z] = (state[x][y] >> (z * 8)) & 0xFF;
            }
        }
    }
}

/*
 * sha3_256: compute SHA3-256 digest of a file.
 *
 * Parameters:
 * - filename: path to input file; opened in binary mode ("rb").
 * - digest: output buffer of 32 bytes.
 *
 * Process:
 * - Absorb RATE-sized chunks with XOR into the first RATE bytes of the state,
 *   applying keccak_f() after each full chunk.
 * - After EOF, pad remaining bytes with 0x06 domain separator and 0x80 final bit,
 *   absorb and apply a final keccak_f().
 * - Squeeze first 32 bytes of the state as the SHA3-256 digest.
 *
 * Returns:
 * - 0 on success; 1 if fopen failed; 2 if fread error occurred.
 */
int sha3_256(const char *filename, uint8_t digest[DIGEST_SIZE]) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        fprintf(stderr, "Error opening file %s\n", filename);
        return 1;
    }

    uint8_t state[STATE_SIZE] = {0};
    uint8_t buffer[RATE];
    size_t len;

    while ((len = fread(buffer, 1, RATE, file)) == RATE) {
        for (size_t i = 0; i < RATE; i++) {
            state[i] ^= buffer[i];
        }
        keccak_f(state);
    }

    if (ferror(file)) {
        fclose(file);
        return 2;
    }

    memset(buffer + len, 0, RATE - len);
    buffer[len] = 0x06;
    buffer[RATE - 1] |= 0x80;
    for (size_t i = 0; i < RATE; i++) {
        state[i] ^= buffer[i];
    }
    keccak_f(state);

    memcpy(digest, state, DIGEST_SIZE);

    fclose(file);
    return 0;
}

/* Print lowercased bytes. */
void print_hex(const uint8_t *data, size_t len) {
    for (size_t i = 0; i < len; i++) {
        printf("%02x", data[i]);
    }
    printf("\n");
}

/*
 * main: compute and print SHA3-256 digest of a file.
 * Usage: <program> <filename>
 * Returns 0 on success, non-zero on error.
 */
int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    uint8_t digest[DIGEST_SIZE];
    int rc = sha3_256(argv[1], digest);
    if (rc != 0) {
        return rc;
    }
    print_hex(digest, DIGEST_SIZE);
    return 0;
}