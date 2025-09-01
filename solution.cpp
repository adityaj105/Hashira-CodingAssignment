#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <stdexcept>


using namespace std;


long double convertToDecimal(const string &val, int base) {
    long double result = 0;
    for (char c : val) {
        int digit;
        if (isdigit(c)) digit = c - '0';
        else digit = tolower(c) - 'a' + 10;
        if (digit >= base) {
            cerr << "Invalid digit for base " << base << ": " << c << "\n";
            exit(1);
        }
        result = result * base + digit;
    }
    return result;
}


vector<long double> lagrangeInterpolation(const vector<long double>& x, const vector<long double>& y) {
    int n = x.size();
    vector<long double> coef(n, 0.0);

    for (int i = 0; i < n; i++) {
        vector<long double> basis = {1.0};
        long double denom = 1.0;

        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            vector<long double> next(basis.size() + 1, 0.0);
            for (int k = 0; k < basis.size(); k++) {
                next[k] -= basis[k] * x[j];
                next[k+1] += basis[k];
            }
            basis = next;
            denom *= (x[i] - x[j]);
        }

        for (int k = 0; k < basis.size(); k++)
            coef[k] += y[i] * basis[k] / denom;
    }

    return coef;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string input, json;
    while (getline(cin, input)) json += input;
    if (json.empty()) {
        cerr << "Empty input\n";
        return 1;
    }

    // Remove spaces
    json.erase(remove(json.begin(), json.end(), ' '), json.end());

    // Extract n and k
    size_t pos = json.find("\"n\":");
    if (pos == string::npos) { cerr << "Missing n\n"; return 1; }
    int n = stoi(json.substr(pos + 4, json.find(',', pos) - (pos + 4)));

    pos = json.find("\"k\":");
    if (pos == string::npos) { cerr << "Missing k\n"; return 1; }
    int k = stoi(json.substr(pos + 4, json.find('}', pos) - (pos + 4)));

    // Collect roots
    vector<long double> X, Y;
    for (int i = 1; i <= n; i++) {
        string key = "\"" + to_string(i) + "\":";
        pos = json.find(key);
        if (pos == string::npos) continue;

        size_t bpos = json.find("\"base\":\"", pos);
        if (bpos == string::npos) continue;
        bpos += 8;
        size_t bend = json.find("\"", bpos + 1);
        int base = stoi(json.substr(bpos + 1, bend - (bpos + 1)));
        if (base < 2 || base > 16) {
            cerr << "Unsupported base: " << base << "\n";
            return 1;
        }

        size_t vpos = json.find("\"value\":\"", pos);
        if (vpos == string::npos) continue;
        vpos += 9;
        size_t vend = json.find("\"", vpos + 1);
        string value = json.substr(vpos + 1, vend - (vpos + 1));

        long double val = convertToDecimal(value, base);
        X.push_back(i);
        Y.push_back(val);
    }

    if (X.size() < k) {
        cerr << "Not enough roots to reconstruct polynomial\n";
        return 1;
    }

    // Use first k roots
    vector<long double> coef = lagrangeInterpolation(
        vector<long double>(X.begin(), X.begin() + k),
        vector<long double>(Y.begin(), Y.begin() + k)
    );

    // Output decimal roots
    cout << "Converted Roots:\n";
    for (auto v : Y) cout << fixed << setprecision(0) << v << " ";
    cout << "\n\n";

    // Output coefficients
    cout << "Polynomial Coefficients (from constant term to degree m):\n";
    for (auto c : coef) cout << fixed << setprecision(2) << c << " ";
    cout << "\n";

    return 0;
}
