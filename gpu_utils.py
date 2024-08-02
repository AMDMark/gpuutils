import cupy as cp
from Crypto.Hash import RIPEMD160
import hashlib
import base58

sha256_kernel_code = '''
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;

__constant__ uint32_t k[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

__device__ uint32_t rotr(uint32_t x, uint32_t n) {
    return (x >> n) | (x << (32 - n));
}

extern "C" __global__ void sha256_transform(const uint8_t *data, uint8_t *hash, uint32_t *w_out) {
    uint32_t h0 = 0x6a09e667;
    uint32_t h1 = 0xbb67ae85;
    uint32_t h2 = 0x3c6ef372;
    uint32_t h3 = 0xa54ff53a;
    uint32_t h4 = 0x510e527f;
    uint32_t h5 = 0x9b05688c;
    uint32_t h6 = 0x1f83d9ab;
    uint32_t h7 = 0x5be0cd19;

    uint32_t w[64];
    for (int i = 0; i < 16; ++i) {
        w[i] = (data[i * 4] << 24) | (data[i * 4 + 1] << 16) | (data[i * 4 + 2] << 8) | (data[i * 4 + 3]);
    }
    for (int i = 16; i < 64; ++i) {
        uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
        uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }

    for (int i = 0; i < 64; ++i) {
        w_out[i] = w[i];
    }

    uint32_t a = h0;
    uint32_t b = h1;
    uint32_t c = h2;
    uint32_t d = h3;
    uint32_t e = h4;
    uint32_t f = h5;
    uint32_t g = h6;
    uint32_t h = h7;

    for (int i = 0; i < 64; ++i) {
        uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        uint32_t ch = (e & f) ^ (~e & g);
        uint32_t temp1 = h + S1 + ch + k[i] + w[i];
        uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t temp2 = S0 + maj;

        h = g;
        g = f;
        f = e;
        e = d + temp1;
        d = c;
        c = b;
        b = a;
        a = temp1 + temp2;
    }

    h0 += a;
    h1 += b;
    h2 += c;
    h3 += d;
    h4 += e;
    h5 += f;
    h6 += g;
    h7 += h;

    hash[0] = (h0 >> 24) & 0xff;
    hash[1] = (h0 >> 16) & 0xff;
    hash[2] = (h0 >> 8) & 0xff;
    hash[3] = h0 & 0xff;
    hash[4] = (h1 >> 24) & 0xff;
    hash[5] = (h1 >> 16) & 0xff;
    hash[6] = (h1 >> 8) & 0xff;
    hash[7] = h1 & 0xff;
    hash[8] = (h2 >> 24) & 0xff;
    hash[9] = (h2 >> 16) & 0xff;
    hash[10] = (h2 >> 8) & 0xff;
    hash[11] = h2 & 0xff;
    hash[12] = (h3 >> 24) & 0xff;
    hash[13] = (h3 >> 16) & 0xff;
    hash[14] = (h3 >> 8) & 0xff;
    hash[15] = h3 & 0xff;
    hash[16] = (h4 >> 24) & 0xff;
    hash[17] = (h4 >> 16) & 0xff;
    hash[18] = (h4 >> 8) & 0xff;
    hash[19] = h4 & 0xff;
    hash[20] = (h5 >> 24) & 0xff;
    hash[21] = (h5 >> 16) & 0xff;
    hash[22] = (h5 >> 8) & 0xff;
    hash[23] = h5 & 0xff;
    hash[24] = (h6 >> 24) & 0xff;
    hash[25] = (h6 >> 16) & 0xff;
    hash[26] = (h6 >> 8) & 0xff;
    hash[27] = h6 & 0xff;
    hash[28] = (h7 >> 24) & 0xff;
    hash[29] = (h7 >> 16) & 0xff;
    hash[30] = (h7 >> 8) & 0xff;
    hash[31] = h7 & 0xff;
}
'''

sha256_kernel = cp.RawKernel(sha256_kernel_code, 'sha256_transform')

def pad_data(data):
    original_length = len(data)
    padding = b'\x80' + b'\x00' * ((56 - (original_length + 1) % 64) % 64)
    padded_data = data + padding + (original_length * 8).to_bytes(8, byteorder='big')
    return padded_data

def sha256_gpu(data):
    padded_data = pad_data(data)
    data_gpu = cp.array(list(padded_data), dtype=cp.uint8)  # Convert bytes to list of integers
    hash_gpu = cp.empty(32, dtype=cp.uint8)
    w_gpu = cp.empty(64, dtype=cp.uint32)

    sha256_kernel((1,), (1,), (data_gpu, hash_gpu, w_gpu))

    hash_result = cp.asnumpy(hash_gpu)
    w_result = cp.asnumpy(w_gpu)
    return hash_result.tobytes(), w_result

def ripemd160_cpu(data):
    ripemd160 = RIPEMD160.new()
    ripemd160.update(data)
    return ripemd160.digest()

def base58check_encode(data):
    checksum = hashlib.sha256(hashlib.sha256(data).digest()).digest()[:4]
    data_with_checksum = data + checksum
    encoded = base58.b58encode(data_with_checksum).decode()
    return encoded
