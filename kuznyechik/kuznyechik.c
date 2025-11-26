#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <windows.h>

#ifndef NTSTATUS
typedef LONG NTSTATUS;
#endif

#ifndef ULONG
typedef unsigned long ULONG;
#endif
#ifndef BCRYPT_USE_SYSTEM_PREFERRED_RNG
#define BCRYPT_USE_SYSTEM_PREFERRED_RNG 0x00000002u
#endif
#if defined(__MINGW32__) || defined(__MINGW64__)
typedef unsigned char UCHAR;
__declspec(dllimport) NTSTATUS WINAPI BCryptGenRandom(void* hAlgorithm, UCHAR* pbBuffer, ULONG cbBuffer, ULONG dwFlags);
#else
#include <bcrypt.h>
#endif

#define KUZNECHIK_BLOCK_SIZE 16

typedef uint64_t chunk[2];

#define ERR_OK 0
#define ERR_KEY_OPEN 1
#define ERR_KEY_SIZE 2
#define ERR_SAVE_OPEN 3
#define ERR_SAVE_WRITE 4
#define ERR_RNG 5
#define ERR_OPEN_IN 6
#define ERR_OPEN_OUT 7
#define ERR_RNG_IV 8
#define ERR_WRITE_IV 9
#define ERR_WRITE_BODY 10
#define ERR_READ_IV 11
#define ERR_HEX_OPEN 12
#define ERR_BAD_ARG 100
#define ERR_PATH_LEN 101

static int validate_path(const char *p) {
    if (!p) {
        return ERR_PATH_LEN;
    }
    size_t n = strlen(p);
    if (n == 0 || n >= MAX_PATH) {
        return ERR_PATH_LEN;
    }
    return ERR_OK;
}

const uint8_t Pi[] = {
    0xFC, 0xEE, 0xDD, 0x11, 0xCF, 0x6E, 0x31, 0x16, 0xFB, 0xC4, 0xFA, 0xDA, 0x23, 0xC5, 0x04, 0x4D,
    0xE9, 0x77, 0xF0, 0xDB, 0x93, 0x2E, 0x99, 0xBA, 0x17, 0x36, 0xF1, 0xBB, 0x14, 0xCD, 0x5F, 0xC1,
    0xF9, 0x18, 0x65, 0x5A, 0xE2, 0x5C, 0xEF, 0x21, 0x81, 0x1C, 0x3C, 0x42, 0x8B, 0x01, 0x8E, 0x4F,
    0x05, 0x84, 0x02, 0xAE, 0xE3, 0x6A, 0x8F, 0xA0, 0x06, 0x0B, 0xED, 0x98, 0x7F, 0xD4, 0xD3, 0x1F,
    0xEB, 0x34, 0x2C, 0x51, 0xEA, 0xC8, 0x48, 0xAB, 0xF2, 0x2A, 0x68, 0xA2, 0xFD, 0x3A, 0xCE, 0xCC,
    0xB5, 0x70, 0x0E, 0x56, 0x08, 0x0C, 0x76, 0x12, 0xBF, 0x72, 0x13, 0x47, 0x9C, 0xB7, 0x5D, 0x87,
    0x15, 0xA1, 0x96, 0x29, 0x10, 0x7B, 0x9A, 0xC7, 0xF3, 0x91, 0x78, 0x6F, 0x9D, 0x9E, 0xB2, 0xB1,
    0x32, 0x75, 0x19, 0x3D, 0xFF, 0x35, 0x8A, 0x7E, 0x6D, 0x54, 0xC6, 0x80, 0xC3, 0xBD, 0x0D, 0x57,
    0xDF, 0xF5, 0x24, 0xA9, 0x3E, 0xA8, 0x43, 0xC9, 0xD7, 0x79, 0xD6, 0xF6, 0x7C, 0x22, 0xB9, 0x03,
    0xE0, 0x0F, 0xEC, 0xDE, 0x7A, 0x94, 0xB0, 0xBC, 0xDC, 0xE8, 0x28, 0x50, 0x4E, 0x33, 0x0A, 0x4A,
    0xA7, 0x97, 0x60, 0x73, 0x1E, 0x00, 0x62, 0x44, 0x1A, 0xB8, 0x38, 0x82, 0x64, 0x9F, 0x26, 0x41,
    0xAD, 0x45, 0x46, 0x92, 0x27, 0x5E, 0x55, 0x2F, 0x8C, 0xA3, 0xA5, 0x7D, 0x69, 0xD5, 0x95, 0x3B,
    0x07, 0x58, 0xB3, 0x40, 0x86, 0xAC, 0x1D, 0xF7, 0x30, 0x37, 0x6B, 0xE4, 0x88, 0xD9, 0xE7, 0x89,
    0xE1, 0x1B, 0x83, 0x49, 0x4C, 0x3F, 0xF8, 0xFE, 0x8D, 0x53, 0xAA, 0x90, 0xCA, 0xD8, 0x85, 0x61,
    0x20, 0x71, 0x67, 0xA4, 0x2D, 0x2B, 0x09, 0x5B, 0xCB, 0x9B, 0x25, 0xD0, 0xBE, 0xE5, 0x6C, 0x52,
    0x59, 0xA6, 0x74, 0xD2, 0xE6, 0xF4, 0xB4, 0xC0, 0xD1, 0x66, 0xAF, 0xC2, 0x39, 0x4B, 0x63, 0xB6
};

const uint8_t Pi_reverse[] = {
    0xA5, 0x2D, 0x32, 0x8F, 0x0E, 0x30, 0x38, 0xC0, 0x54, 0xE6, 0x9E, 0x39, 0x55, 0x7E, 0x52, 0x91,
    0x64, 0x03, 0x57, 0x5A, 0x1C, 0x60, 0x07, 0x18, 0x21, 0x72, 0xA8, 0xD1, 0x29, 0xC6, 0xA4, 0x3F,
    0xE0, 0x27, 0x8D, 0x0C, 0x82, 0xEA, 0xAE, 0xB4, 0x9A, 0x63, 0x49, 0xE5, 0x42, 0xE4, 0x15, 0xB7,
    0xC8, 0x06, 0x70, 0x9D, 0x41, 0x75, 0x19, 0xC9, 0xAA, 0xFC, 0x4D, 0xBF, 0x2A, 0x73, 0x84, 0xD5,
    0xC3, 0xAF, 0x2B, 0x86, 0xA7, 0xB1, 0xB2, 0x5B, 0x46, 0xD3, 0x9F, 0xFD, 0xD4, 0x0F, 0x9C, 0x2F,
    0x9B, 0x43, 0xEF, 0xD9, 0x79, 0xB6, 0x53, 0x7F, 0xC1, 0xF0, 0x23, 0xE7, 0x25, 0x5E, 0xB5, 0x1E,
    0xA2, 0xDF, 0xA6, 0xFE, 0xAC, 0x22, 0xF9, 0xE2, 0x4A, 0xBC, 0x35, 0xCA, 0xEE, 0x78, 0x05, 0x6B,
    0x51, 0xE1, 0x59, 0xA3, 0xF2, 0x71, 0x56, 0x11, 0x6A, 0x89, 0x94, 0x65, 0x8C, 0xBB, 0x77, 0x3C,
    0x7B, 0x28, 0xAB, 0xD2, 0x31, 0xDE, 0xC4, 0x5F, 0xCC, 0xCF, 0x76, 0x2C, 0xB8, 0xD8, 0x2E, 0x36,
    0xDB, 0x69, 0xB3, 0x14, 0x95, 0xBE, 0x62, 0xA1, 0x3B, 0x16, 0x66, 0xE9, 0x5C, 0x6C, 0x6D, 0xAD,
    0x37, 0x61, 0x4B, 0xB9, 0xE3, 0xBA, 0xF1, 0xA0, 0x85, 0x83, 0xDA, 0x47, 0xC5, 0xB0, 0x33, 0xFA,
    0x96, 0x6F, 0x6E, 0xC2, 0xF6, 0x50, 0xFF, 0x5D, 0xA9, 0x8E, 0x17, 0x1B, 0x97, 0x7D, 0xEC, 0x58,
    0xF7, 0x1F, 0xFB, 0x7C, 0x09, 0x0D, 0x7A, 0x67, 0x45, 0x87, 0xDC, 0xE8, 0x4F, 0x1D, 0x4E, 0x04,
    0xEB, 0xF8, 0xF3, 0x3E, 0x3D, 0xBD, 0x8A, 0x88, 0xDD, 0xCD, 0x0B, 0x13, 0x98, 0x02, 0x93, 0x80,
    0x90, 0xD0, 0x24, 0x34, 0xCB, 0xED, 0xF4, 0xCE, 0x99, 0x10, 0x44, 0x40, 0x92, 0x3A, 0x01, 0x26,
    0x12, 0x1A, 0x48, 0x68, 0xF5, 0x81, 0x8B, 0xC7, 0xD6, 0x20, 0x0A, 0x08, 0x00, 0x4C, 0xD7, 0x74
};

const uint8_t linear_vector[] = {
    0x94, 0x20, 0x85, 0x10, 0xC2, 0xC0, 0x01, 0xFB,
    0x01, 0xC0, 0xC2, 0x10, 0x85, 0x20, 0x94, 0x01
};

void X(chunk a, chunk b, chunk c) {
    c[0] = a[0] ^ b[0];
    c[1] = a[1] ^ b[1];
}

void S(chunk in_out) {
    int i;
    uint8_t *byte = (int8_t*) in_out;
    for (i = 0; i < KUZNECHIK_BLOCK_SIZE; i++) {
        byte[i] = Pi[byte[i]];
    }
}

void S_reverse(chunk in_out) {
    int i;
    uint8_t *byte = (int8_t*) in_out;
    for (i = 0; i < KUZNECHIK_BLOCK_SIZE; i++) {
        byte[i] = Pi_reverse[byte[i]];
    }
}

uint8_t GF_mult(uint8_t a, uint8_t b) {    
    uint8_t c;
    
    c = 0;
    while (b) {        
        if (b & 1)
            c ^= a;
        a = (a << 1) ^ (a & 0x80 ? 0xC3 : 0x00);
        b >>= 1;
    }
        
    return c;
}

void R(uint8_t *in_out) {
    int i;
    uint8_t acc = in_out[15];
    uint8_t *byte = (int8_t*) in_out;
    for (i = 14; i >= 0; i--) {
        byte[i + 1] = byte[i];
        acc ^= GF_mult(byte[i], linear_vector[i]);
    }
    byte[0] = acc;
}

void R_reverse(uint8_t *in_out) {
    int i;
    uint8_t acc = in_out[0];
    uint8_t *byte = (int8_t*) in_out;

    for (int i = 0; i < 15; i++) {
        byte[i] = byte[i + 1];
        acc ^= GF_mult(byte[i], linear_vector[i]);
    }

    byte[15] = acc;
}

void L(uint8_t* in_out) {
    int i;
    for (i = 0; i < KUZNECHIK_BLOCK_SIZE; i++) {
        R(in_out);
    }
}

void L_reverse(uint8_t *in_out) {
    int i;
    for (int i = 0; i < 16; i++) {
        R_reverse(in_out);
    }
}

void gen_round_keys(uint8_t* key, chunk* round_keys) {
    int i;

    uint8_t cs[32][KUZNECHIK_BLOCK_SIZE] = {};

    for (i = 0; i < 32; i++) {
        cs[i][15] = i + 1;
        L(cs[i]);
    }

    chunk ks[2] = {};
    round_keys[0][0] = ks[0][0] = ((chunk*) key)[0][0];
    round_keys[0][1] = ks[0][1] = ((chunk*) key)[0][1];
    round_keys[1][0] = ks[1][0] = ((chunk*) key)[1][0];
    round_keys[1][1] = ks[1][1] = ((chunk*) key)[1][1];

    for (i = 1; i <= 32; i++) {
        chunk new_key = {0};

        X(ks[0], (void*)cs[i - 1], new_key);
        S(new_key);
        L((uint8_t*)&new_key);
        X(new_key, ks[1], new_key);

        ks[1][0] = ks[0][0];
        ks[1][1] = ks[0][1];

        ks[0][0] = new_key[0];
        ks[0][1] = new_key[1];

        if ((i > 0) && (i % 8 == 0)) {
            round_keys[(i >> 2)][0] = ks[0][0];
            round_keys[(i >> 2)][1] = ks[0][1];

            round_keys[(i >> 2) + 1][0] = ks[1][0];
            round_keys[(i >> 2) + 1][1] = ks[1][1];
        }
    }
}

void kuznechik_encrypt(chunk *round_keys, chunk in, chunk out) {
    int i;
    chunk p;
    memcpy(p, in, sizeof(chunk));
    for (i = 0; i < 10; i++) {
        X(p, round_keys[i], p);
        if (i < 9) {
            S(p);
            L((uint8_t*)&p);
        }
    }
    memcpy(out, p, sizeof(chunk));
}

void kuznechik_decrypt(chunk *round_keys, chunk in, chunk out) {
    int i;
    chunk p;
    memcpy(p, in, sizeof(chunk));

    X(p, round_keys[9], p);
    for (i = 8; i >= 0; i--) {
        L_reverse((uint8_t*)&p);
        S_reverse(p);
        X(p, round_keys[i], p);
    }
    memcpy(out, p, sizeof(chunk));
}

void print_chunk(chunk p) {
    int i;
    for (i = 0; i < sizeof(chunk); i++) {
        printf("0x%02X ", ((uint8_t *)p)[i]);
    }    
    printf("\n");
}

static void incr128(uint8_t *ctr) {
    for (int i = 15; i >= 0; --i) {
        if (++ctr[i] != 0) break;
    }
}

static void block_encrypt_bytes(chunk *round_keys, const uint8_t in[16], uint8_t out[16]) {
    chunk blk_in;
    chunk blk_out;
    memcpy(&blk_in, in, 16);
    kuznechik_encrypt(round_keys, blk_in, blk_out);
    memcpy(out, &blk_out, 16);
}

static int load_key(const char *path, uint8_t key[32]) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        return ERR_KEY_OPEN;
    }
    size_t n = fread(key, 1, 32, f);
    fclose(f);
    if (n != 32) {
        return ERR_KEY_SIZE;
    }
    return ERR_OK;
}

static int save_key(const char *path, const uint8_t key[32]) {
    FILE *f = fopen(path, "wb");
    if (!f) {
        return ERR_SAVE_OPEN;
    }
    size_t n = fwrite(key, 1, 32, f);
    fclose(f);
    if (n != 32) {
        return ERR_SAVE_WRITE;
    }
    return ERR_OK;
}

static int genkey_cmd(const char *path) {
    uint8_t key[32];
    NTSTATUS st = BCryptGenRandom(NULL, key, (ULONG)sizeof(key), BCRYPT_USE_SYSTEM_PREFERRED_RNG);
    if (st != 0) {
        return ERR_RNG;
    }
    return save_key(path, key);
}

static int enc_ctr(const char *key_path, const char *in_path, const char *out_path) {
    uint8_t key[32];
    int lk = load_key(key_path, key);
    if (lk != ERR_OK) {
        return lk;
    }

    chunk round_keys[10] = {};
    gen_round_keys(key, round_keys);

    FILE *fi = fopen(in_path, "rb");
    if (!fi) {
        return ERR_OPEN_IN;
    }
    FILE *fo = fopen(out_path, "wb");
    if (!fo) {
        fclose(fi);
        return ERR_OPEN_OUT;
    }

    uint8_t iv[16];
    if (BCryptGenRandom(NULL, iv, sizeof iv, BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0) {
        fclose(fi);
        fclose(fo);
        return ERR_RNG_IV;
    }
    if (fwrite(iv, 1, 16, fo) != 16) {
        fclose(fi);
        fclose(fo);
        return ERR_WRITE_IV;
    }

    uint8_t ctr[16];
    memcpy(ctr, iv, 16);

    uint8_t inbuf[65536];
    uint8_t keystream[16];
    size_t r;
    do {
        r = fread(inbuf, 1, sizeof inbuf, fi);
        if (r > 0) {
            size_t off = 0;
            while (off < r) {
                block_encrypt_bytes(round_keys, ctr, keystream);
                size_t n = (r - off) >= 16 ? 16 : (r - off);
                for (size_t i = 0; i < n; ++i) inbuf[off + i] ^= keystream[i];
                incr128(ctr);
                off += n;
            }
            if (fwrite(inbuf, 1, r, fo) != r) {
                fclose(fi);
                fclose(fo);
                return ERR_WRITE_BODY;
            }
        }
    } while (r > 0);

    fclose(fi);
    fclose(fo);
    return 0;
}

static int dec_ctr(const char *key_path, const char *in_path, const char *out_path) {
    uint8_t key[32];
    if (load_key(key_path, key) != 0) {
        return 1;
    }

    chunk round_keys[10] = {};
    gen_round_keys(key, round_keys);

    FILE *fi = fopen(in_path, "rb");
    if (!fi) {
        return 2;
    }
    FILE *fo = fopen(out_path, "wb");
    if (!fo) {
        fclose(fi);
        return 3;
    }

    uint8_t iv[16];
    if (fread(iv, 1, 16, fi) != 16) {
        fclose(fi);
        fclose(fo);
        return ERR_READ_IV;
    }

    uint8_t ctr[16];
    memcpy(ctr, iv, 16);

    uint8_t inbuf[65536];
    uint8_t keystream[16];
    size_t r;
    do {
        r = fread(inbuf, 1, sizeof inbuf, fi);
        if (r > 0) {
            size_t off = 0;
            while (off < r) {
                block_encrypt_bytes(round_keys, ctr, keystream);
                size_t n = (r - off) >= 16 ? 16 : (r - off);
                for (size_t i = 0; i < n; ++i) inbuf[off + i] ^= keystream[i];
                incr128(ctr);
                off += n;
            }
            if (fwrite(inbuf, 1, r, fo) != r) {
                fclose(fi);
                fclose(fo);
                return ERR_WRITE_BODY;
            }
        }
    } while (r > 0);

    fclose(fi);
    fclose(fo);
    return 0;
}

static int hex_dump(const char *path)
{
    FILE *f = fopen(path, "rb");
    if (!f) {
        return ERR_HEX_OPEN;
    }
    int c;
    size_t i = 0;
    while ((c = fgetc(f)) != EOF) {
        printf("%02X", (unsigned)(uint8_t)c);
        if (++i % 1 == 0) printf(" ");
    }

    printf("\n");
    fclose(f);
    return 0;
}

static int file_exists_regular(const char *p)
{
    DWORD a = GetFileAttributesA(p);
    if (a == INVALID_FILE_ATTRIBUTES)
    {
        return 0;
    }
    if (a & FILE_ATTRIBUTE_DIRECTORY)
    {
        return 0;
    }
    return 1;
}

static void print_usage(const char *prog)
{
    fprintf(stderr,
        "Usage:\n"
        "  %s genkey <key.bin>\n"
        "  %s enc <key.bin> <in> <out>\n"
        "  %s dec <key.bin> <in> <out>\n"
        "  %s hex <file>\n",
        prog, prog, prog, prog);
}

static void explain_err(int code)
{
    switch (code)
    {
        case ERR_KEY_OPEN: fprintf(stderr, "error: cannot open key file\n"); break;
        case ERR_KEY_SIZE: fprintf(stderr, "error: key file must be exactly 32 bytes\n"); break;
        case ERR_OPEN_IN: fprintf(stderr, "error: cannot open input file\n"); break;
        case ERR_OPEN_OUT: fprintf(stderr, "error: cannot open output file\n"); break;
        case ERR_RNG: fprintf(stderr, "error: RNG failure\n"); break;
        case ERR_RNG_IV: fprintf(stderr, "error: cannot generate IV\n"); break;
        case ERR_READ_IV: fprintf(stderr, "error: cannot read IV from input\n"); break;
        case ERR_WRITE_IV: fprintf(stderr, "error: cannot write IV to output\n"); break;
        case ERR_WRITE_BODY: fprintf(stderr, "error: write failed\n"); break;
        case ERR_HEX_OPEN: fprintf(stderr, "error: cannot open file for hex dump\n"); break;
        case ERR_PATH_LEN: fprintf(stderr, "error: bad path length\n"); break;
        case ERR_BAD_ARG: fprintf(stderr, "error: bad arguments\n"); break;
        default: break;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage(argv[0]);
        return ERR_BAD_ARG;
    }
    if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        print_usage(argv[0]);
        return 0;
    }
    if (strcmp(argv[1], "genkey") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Usage: %s genkey <key.bin>\n", argv[0]);
            return ERR_BAD_ARG;
        }
        int vp = validate_path(argv[2]);
        if (vp) {
            fprintf(stderr, "bad path length: %s\n", argv[2]);
            return vp;
        }
        int rc = genkey_cmd(argv[2]);
        if (rc) {
            fprintf(stderr, "genkey failed (%d)\n", rc);
            explain_err(rc);
            return rc;
        }
        return 0;
    }
    else if (strcmp(argv[1], "enc") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Usage: %s enc <key.bin> <in> <out>\n", argv[0]);
            return ERR_BAD_ARG;
        }
        int v1 = validate_path(argv[2]);
        int v2 = validate_path(argv[3]);
        int v3 = validate_path(argv[4]);
        if (v1||v2||v3) {
            fprintf(stderr, "bad path length\n");
            return ERR_PATH_LEN;
        }
        if (strcmp(argv[3], argv[4]) == 0) {
            fprintf(stderr, "error: input and output paths must differ\n");
            return ERR_BAD_ARG;
        }
        if (strcmp(argv[2], argv[4]) == 0) {
            fprintf(stderr, "error: key path must differ from output path\n");
            return ERR_BAD_ARG;
        }
        if (!file_exists_regular(argv[3])) {
            fprintf(stderr, "error: input file does not exist or is a directory\n");
            return ERR_OPEN_IN;
        }
        int rc = enc_ctr(argv[2], argv[3], argv[4]);
        if (rc) {
            fprintf(stderr, "enc failed (%d)\n", rc);
            explain_err(rc);
            return rc;
        }
        return 0;
    }
    else if (strcmp(argv[1], "dec") == 0) {
        if (argc != 5) {
            fprintf(stderr, "Usage: %s dec <key.bin> <in> <out>\n", argv[0]);
            return ERR_BAD_ARG;
        }
        int v1 = validate_path(argv[2]);
        int v2 = validate_path(argv[3]);
        int v3 = validate_path(argv[4]);
        if (v1||v2||v3) {
            fprintf(stderr, "bad path length\n");
            return ERR_PATH_LEN;
        }
        if (strcmp(argv[3], argv[4]) == 0) {
            fprintf(stderr, "error: input and output paths must differ\n");
            return ERR_BAD_ARG;
        }
        if (strcmp(argv[2], argv[4]) == 0) {
            fprintf(stderr, "error: key path must differ from output path\n");
            return ERR_BAD_ARG;
        }
        if (!file_exists_regular(argv[3])) {
            fprintf(stderr, "error: input file does not exist or is a directory\n");
            return ERR_OPEN_IN;
        }
        int rc = dec_ctr(argv[2], argv[3], argv[4]);
        if (rc) {
            fprintf(stderr, "dec failed (%d)\n", rc);
            explain_err(rc);
            return rc;
        }
        return 0;
    }
    else if (strcmp(argv[1], "hex") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Usage: %s hex <file>\n", argv[0]);
            return ERR_BAD_ARG;
        }
        int vp = validate_path(argv[2]);
        if (vp) {
            fprintf(stderr, "bad path length: %s\n", argv[2]);
            return vp;
        }
        if (!file_exists_regular(argv[2])) {
            fprintf(stderr, "error: file does not exist or is a directory\n");
            return ERR_HEX_OPEN;
        }
        int rc = hex_dump(argv[2]);
        if (rc) {
            fprintf(stderr, "hex failed (%d)\n", rc);
            explain_err(rc);
            return rc;
        }
        return 0;
    }
    else {
        fprintf(stderr, "Unknown command: %s\n", argv[1]);
        print_usage(argv[0]);
        return ERR_BAD_ARG;
    }
}