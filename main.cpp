#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <queue>
#include <chrono>
#include <random>
#include <cassert>

using namespace std;

pair<vector<uint8_t>, int> bwtTransform(const vector<uint8_t>& s) {
    int n = s.size();
    if (n == 0) return {{}, 0};

    vector<int> idx(n);
    for (int i = 0; i < n; i++) idx[i] = i;

    sort(idx.begin(), idx.end(), [&](int a, int b) {
        for (int k = 0; k < n; k++) {
            uint8_t ca = s[(a + k) % n];
            uint8_t cb = s[(b + k) % n];
            if (ca != cb) return ca < cb;
        }
        return false;
    });

    vector<uint8_t> last(n);
    int primaryIndex = 0;

    for (int i = 0; i < n; i++) {
        int j = idx[i];
        last[i] = s[(j + n - 1) % n];
        if (j == 0) primaryIndex = i;
    }

    return {last, primaryIndex};
}

vector<uint8_t> inverseBWT(const vector<uint8_t>& L, int primaryIndex) {
    int n = L.size();
    if (n == 0) return {};

    vector<int> count(256, 0);
    for (auto c : L) count[c]++;

    vector<int> start(256, 0);
    for (int i = 1; i < 256; i++)
        start[i] = start[i - 1] + count[i - 1];

    vector<int> T(n);
    vector<int> occ(256, 0);

    for (int i = 0; i < n; i++) {
        uint8_t c = L[i];
        T[i] = start[c] + occ[c]++;
    }

    vector<uint8_t> result(n);
    int pos = primaryIndex;

    for (int i = n - 1; i >= 0; i--) {
        result[i] = L[pos];
        pos = T[pos];
    }

    return result;
}


vector<int> mtfEncode(const vector<uint8_t>& data) {
    vector<int> res;
    vector<uint8_t> alphabet(256);

    for (int i = 0; i < 256; i++) alphabet[i] = i;

    for (auto c : data) {
        int pos = 0;
        while (alphabet[pos] != c) pos++;

        res.push_back(pos);

        alphabet.erase(alphabet.begin() + pos);
        alphabet.insert(alphabet.begin(), c);
    }

    return res;
}

vector<uint8_t> mtfDecode(const vector<int>& data) {
    vector<uint8_t> res;
    vector<uint8_t> alphabet(256);

    for (int i = 0; i < 256; i++) alphabet[i] = i;

    for (int idx : data) {
        uint8_t c = alphabet[idx];
        res.push_back(c);

        alphabet.erase(alphabet.begin() + idx);
        alphabet.insert(alphabet.begin(), c);
    }

    return res;
}


vector<pair<uint8_t, int>> rleEncodeBytes(const vector<uint8_t>& data) {
    vector<pair<uint8_t, int>> res;

    for (size_t i = 0; i < data.size();) {
        uint8_t c = data[i];
        size_t j = i;

        while (j < data.size() && data[j] == c) j++;

        res.push_back({c, (int)(j - i)});
        i = j;
    }

    return res;
}


vector<pair<uint8_t, int>> rleEncodeMTF(const vector<int>& data) {
    vector<pair<uint8_t, int>> res;

    for (size_t i = 0; i < data.size();) {
    int val = data[i];
assert(val >= 0 && val < 256);

size_t j = i;

while (j < data.size() && data[j] == val) j++;

res.push_back({(uint8_t)val, (int)(j - i)});
        i = j;
    }

    return res;
}

vector<uint8_t> rleDecode(const vector<pair<uint8_t, int>>& data) {
    vector<uint8_t> res;

    for (auto& p : data) {
        for (int i = 0; i < p.second; i++)
            res.push_back(p.first);
    }

    return res;
}

struct Node {
    uint8_t ch;
    int freq;
    Node* left;
    Node* right;

    Node(uint8_t c, int f) : ch(c), freq(f), left(nullptr), right(nullptr) {}
};

struct Compare {
    bool operator()(Node* a, Node* b) {
        return a->freq > b->freq;
    }
};

void buildCodes(Node* root, string code, map<uint8_t, string>& table) {
    if (!root) return;

    if (!root->left && !root->right)
        table[root->ch] = code;

    buildCodes(root->left, code + "0", table);
    buildCodes(root->right, code + "1", table);
}

void deleteTree(Node* root) {
    if (!root) return;
    deleteTree(root->left);
    deleteTree(root->right);
    delete root;
}

int huffmanCompress(const vector<uint8_t>& data) {
    if (data.empty()) return 0;

    map<uint8_t, int> freq;
    for (auto c : data) freq[c]++;

    priority_queue<Node*, vector<Node*>, Compare> pq;

    for (auto& p : freq)
        pq.push(new Node(p.first, p.second));

    while (pq.size() > 1) {
        Node* a = pq.top(); pq.pop();
        Node* b = pq.top(); pq.pop();

        Node* merged = new Node(0, a->freq + b->freq);
        merged->left = a;
        merged->right = b;

        pq.push(merged);
    }

    Node* root = pq.top();

    map<uint8_t, string> table;
    buildCodes(root, "", table);

    int bits = 0;
    for (auto c : data)
        bits += table[c].size();

    deleteTree(root);

    return bits;
}

vector<int> lzwCompress(const vector<uint8_t>& input) {
    map<vector<uint8_t>, int> dict;

    for (int i = 0; i < 256; i++)
        dict[{(uint8_t)i}] = i;

    vector<uint8_t> w;
    vector<int> result;
    int code = 256;

    for (auto c : input) {
        auto wc = w;
        wc.push_back(c);

        if (dict.count(wc)) {
            w = wc;
        } else {
            result.push_back(dict[w]);
            dict[wc] = code++;
            w = {c};
        }
    }

    if (!w.empty())
        result.push_back(dict[w]);

    return result;
}


struct Options {
    bool bwt = false;
    bool mtf = false;
    bool rle = false;
};

struct EncodedData {
    vector<uint8_t> data;
    vector<int> mtf;
    vector<pair<uint8_t,int>> rle;
    int index = 0;
};

EncodedData encodePipeline(const vector<uint8_t>& input, Options opt) {
    EncodedData out;

    vector<uint8_t> data = input;

    // BWT
    if (opt.bwt) {
        auto res = bwtTransform(data);
        data = res.first;
        out.index = res.second;
    }

    // MTF
    if (opt.mtf) {
        out.mtf = mtfEncode(data);
    }

    // RLE
    if (opt.rle) {
        if (opt.mtf) {
            out.rle = rleEncodeMTF(out.mtf);
        } else {
            out.rle = rleEncodeBytes(data);
        }
    }

    if (!opt.mtf && !opt.rle)
        out.data = data;

    return out;
}

vector<uint8_t> decodePipeline(const EncodedData& enc, Options opt) {
    vector<uint8_t> data;

    if (opt.rle) {
        vector<uint8_t> temp = rleDecode(enc.rle);

        if (opt.mtf) {
            vector<int> mtfData;
            for (auto c : temp)
                mtfData.push_back((int)c);

            data = mtfDecode(mtfData);
        } else {
            data = temp;
        }
    }
    else if (opt.mtf) {
        data = mtfDecode(enc.mtf);
    }
    else {
        data = enc.data;
    }

    if (opt.bwt) {
        data = inverseBWT(data, enc.index);
    }

    return data;
}



vector<uint8_t> makeData(const string& s) {
    return vector<uint8_t>(s.begin(), s.end());
}

vector<uint8_t> randomData(int n) {
    vector<uint8_t> v(n);
    random_device rd;
    for (int i = 0; i < n; i++)
        v[i] = rd() % 256;
    return v;
}

void testCase(const string& name, vector<uint8_t> data, Options opt) {
    auto start = chrono::high_resolution_clock::now();

    auto enc = encodePipeline(data, opt);

    vector<uint8_t> processed;

    if (opt.rle) {
        processed = rleDecode(enc.rle);
    } else if (opt.mtf) {
        processed = mtfDecode(enc.mtf);
    } else {
        processed = enc.data;
    }

    int huff = huffmanCompress(processed);
    int lzw = lzwCompress(processed).size() * sizeof(int);

    auto restored = decodePipeline(enc, opt);

    auto end = chrono::high_resolution_clock::now();

    cout << "\n--- " << name << " ---\n";
    cout << "Correct: " << (restored == data ? "YES" : "NO") << endl;
    cout << "Huffman bits: " << huff << endl;
    cout << "LZW bytes: " << lzw << endl;
    cout << "Time: "
         << chrono::duration_cast<chrono::microseconds>(end - start).count()
         << " us\n";
}

int main() {
    auto text = makeData("banana banana banana");
    auto rep = makeData(string(200, 'A'));
    auto rnd = randomData(200);

    Options none{};
    Options bwt{true,false,false};
    Options bwt_mtf{true,true,false};
    Options full{true,true,true};

    cout << "=== NONE ===\n";
    testCase("text", text, none);
    testCase("repeat", rep, none);
    testCase("random", rnd, none);

    cout << "\n=== BWT ===\n";
    testCase("text", text, bwt);

    cout << "\n=== BWT+MTF ===\n";
    testCase("text", text, bwt_mtf);

    cout << "\n=== BWT+MTF+RLE ===\n";
    testCase("text", text, full);

    return 0;
}

