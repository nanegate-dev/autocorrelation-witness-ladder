// Starting program for the autocorrelation problem of AlphaEvolve's repository
// of problems, number 6 (related to difference bases).
// Contract: candidate <input.json> <output.json> --seed N
// Output: {"heights":[h0,h1,...],"floor":p}
//
// The arithmetic is `problem.txt`'s and is not restated here: read it there.
// What matters to this file is that the score is a MINIMUM over the shifts,
// that it is MAXIMISED, and that the solution states the floor it claims
// rather than having it computed for it.
//
// THE DIRECTION IS THE OPPOSITE OF THE SISTER EXAMPLE'S. Beside this one lives
// `alphaevolve-autocorrelation`, problem 2 of the same repository: there the
// constant is the largest that works, a construction gives an UPPER bound, and
// the search MINIMISES a maximum. Here the constant is the smallest that
// works, a construction gives a LOWER bound, and the search MAXIMISES a
// minimum. The two programs look alike and mean opposite things; a sign copied
// from one into the other is the mistake to expect.
//
// WHY THIS METHOD. The score is a MINIMUM, so moving one height changes
// nothing at all until it changes which shift is the argmin, and then it
// changes by a jump. A search with no gradient wanders such a surface. So the
// method is the one that works on the sister problem, mirrored:
//
//   * SMOOTH THE MINIMUM. min_k A_k is replaced by a softmin-weighted average
//     of the A_k at inverse temperature beta. That has a gradient everywhere,
//     and as beta grows it becomes the true minimum.
//   * ANNEAL beta upward. The flat vector is a stationary point of the
//     smoothed objective, so a search that stops warm sits on it and reports
//     the single-interval score.
//   * MOMENTUM, because the smoothed surface is a long narrow valley.
//   * PROJECT back onto the simplex after every step. The objective is
//     invariant under scale, so fixing the total removes a redundant direction
//     and makes the step size mean something; non-negativity is condition 2
//     and the projection enforces it exactly rather than by clipping.
//   * RESTART, because the basins are real.
//
// WHAT THE NUMBERS ARE. A single interval of length L scores (L-1)/L^2, which
// is 0.25 at its best — every height equal is that, and it is the number to
// beat. Published bounds are 0.37 <= C <= 0.411, and the paper collecting
// these problems records that a numerical attack "appears to be difficult".
//
// THE BUDGET IS ITERATIONS, NOT SECONDS. A wall-clock budget would make this
// program answer differently on a fast machine than a slow one, and ADR 0027's
// cache treats (binary, instance, seed) as one answer. So the anneal runs on
// the iteration index and the knobs below are counts.
//
// The score is never printed or written: the system computes it.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::size_t;
using std::string;
using std::vector;

namespace {

// --- reading the instance --------------------------------------------------

struct Reader {
    string text;
    size_t at = 0;

    bool seek(const string& key) {
        const string needle = "\"" + key + "\"";
        const size_t found = text.find(needle);
        if (found == string::npos) return false;
        at = text.find(':', found + needle.size());
        if (at == string::npos) return false;
        ++at;
        while (at < text.size() && std::isspace(static_cast<unsigned char>(text[at]))) ++at;
        return true;
    }

    double number() {
        char* end = nullptr;
        const double value = std::strtod(text.c_str() + at, &end);
        at = static_cast<size_t>(end - text.c_str());
        return value;
    }
};

struct Instance {
    int pieces = 0;
    int per_unit = 0;
};

Instance read_instance(const char* path) {
    std::ifstream file(path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    Reader reader{buffer.str(), 0};
    Instance in;
    if (reader.seek("pieces")) in.pieces = static_cast<int>(reader.number());
    if (reader.seek("per_unit")) in.per_unit = static_cast<int>(reader.number());
    return in;
}

// --- the score --------------------------------------------------------------

// The self-correlation at every shift 0..lags: c[k] = sum over i of h[i]*h[i+k],
// keeping the i where i+k is inside the list.
//
// CORRELATION, not convolution: the pair at shift k is (i, i + k). The sister
// example pairs (i, k - i) instead, which is a different quantity — that one
// is a convolution and its value at shift zero is not the sum of squares.
void correlate(const vector<double>& h, int lags, vector<double>& c) {
    const int n = static_cast<int>(h.size());
    for (int k = 0; k <= lags; ++k) {
        double total = 0.0;
        for (int i = 0; i + k < n; ++i) total += h[i] * h[i + k];
        c[k] = total;
    }
}

double floor_of(const vector<double>& c, int lags) {
    double worst = c[0];
    for (int k = 1; k <= lags; ++k) worst = std::min(worst, c[k]);
    return worst;
}

// per_unit * floor / (sum h)^2, the objective as `problem.txt` states it.
double ratio(const vector<double>& h, double claimed, int per_unit) {
    double mass = 0.0;
    for (const double value : h) mass += value;
    if (mass <= 0.0) return 0.0;
    return per_unit * claimed / (mass * mass);
}

// --- projection onto the simplex --------------------------------------------

// The nearest point with non-negative entries summing to one (Duchi et al.).
// Sort descending, walk the prefix sums for the threshold that makes the
// remainder sum to one, subtract it, clamp at zero.
void project(vector<double>& h, vector<double>& sorted) {
    const int n = static_cast<int>(h.size());
    sorted = h;
    std::sort(sorted.begin(), sorted.end(), std::greater<double>());
    double running = 0.0;
    double theta = 0.0;
    for (int i = 0; i < n; ++i) {
        running += sorted[i];
        const double candidate = (running - 1.0) / (i + 1);
        if (sorted[i] - candidate > 0.0) theta = candidate;
    }
    for (int i = 0; i < n; ++i) h[i] = std::max(0.0, h[i] - theta);
}

// --- randomness -------------------------------------------------------------

struct Rng {
    std::uint64_t state;
    explicit Rng(std::uint64_t seed)
        : state(seed * 6364136223846793005ULL + 1442695040888963407ULL) {}
    std::uint64_t next() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return state;
    }
    double unit() { return static_cast<double>(next() >> 11) * (1.0 / 9007199254740992.0); }
};

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) return 2;
    std::uint64_t seed = 1;
    for (int i = 3; i + 1 < argc; ++i) {
        if (string(argv[i]) == "--seed") seed = std::strtoull(argv[i + 1], nullptr, 10);
    }
    const Instance in = read_instance(argv[1]);
    if (in.pieces <= 0 || in.per_unit <= 0) return 3;

    // THE TWO COUNTS BELOW MULTIPLY, AND THEIR TOP CORNER IS THE DEADLINE.
    // §15 gives a candidate 120 seconds. Whoever widens one of these ranges is
    // widening a product, so the corner was measured rather than assumed — at
    // the shipped instance, on 2026-09-07:
    //
    //     restarts  steps    seconds   score
    //           24  60000        8.5   0.403970   <- the defaults
    //           24  100000      14.6   0.403291
    //           24  200000      29.5   0.403992
    //           64  60000       23.1   0.403970
    //           64  200000      79.0   (the old corner)
    //
    // NEITHER COUNT BUYS ANYTHING, and that is why `steps` stops at 100000
    // rather than at the 200000 it was first written with. Tripling the steps
    // moves the score by 2e-5 and nearly tripling the restarts moves it by
    // nothing at all, while the old corner cost 79 of the 120 seconds. It never
    // timed out — it was worse than that. The objective does not price time, so
    // a variant that drifts to the corner is promoted for a rounding error and
    // then returns five candidates per branch where the defaults return fifty.
    // The ceiling is set by what has failed to pay, not by what fits.
    //
    // The measured corner is now 64 x 100000, about 38 seconds.
    //
    // One step costs about (lags+1) * pieces multiply-adds, under 3200 here.
    constexpr int kRestarts = 61;  // evolve:param restarts int 1 64
    // Gradient steps per start. The anneal is spread over exactly this many,
    // so a shorter run is a faster anneal and not a truncated one.
    constexpr int kSteps = 90255;  // evolve:param steps int 200 100000
    // Step size, before the 1/(1 + decay*t) taper below.
    constexpr double kRate = 0.0001;  // evolve:param rate real 0.0001 1.0
    constexpr double kRateSpread = 8.541015625;  // evolve:param rate_spread real 1.0 100.0
    // How much of the previous step is carried. The smoothed surface is a long
    // narrow valley and this is what walks along it rather than across it.
    constexpr double kMomentum = 0.9882713671875001;  // evolve:param momentum real 0.0 0.99999
    // Where the anneal starts. Low is smooth and leads somewhere.
    constexpr double kBetaFrom = 13868.048828125;  // evolve:param beta_from real 1.0 100000.0
    // Where it ends. This one must be large: the flat vector is a stationary
    // point of the smoothed objective, so an anneal that stops warm sits on it
    // and reports the single-interval score.
    constexpr double kBetaTo = 76855700.1953125;  // evolve:param beta_to real 1000.0 100000000.0

    const int n = in.pieces;
    const int lags = in.per_unit;  // `problem.txt`: the lags run 0..per_unit
    Rng rng(seed);

    vector<double> h(n), c(lags + 1), soft(lags + 1), gradient(n), velocity(n), sorted(n);
    vector<double> best(n, 1.0 / n);
    correlate(best, lags, c);
    double best_score = ratio(best, floor_of(c, lags), in.per_unit);

    for (int restart = 0; restart < kRestarts; ++restart) {
        // The first start is flat, so the search can never do worse than the
        // single interval everyone already knows about; the rest are random.
        if (restart == 0) {
            std::fill(h.begin(), h.end(), 1.0 / n);
        } else {
            // Start from the best solution found so far, perturbed along a random direction.
            double mass = 0.0;
            double norm = 0.0;
            for (int i = 0; i < n; ++i) { double r = rng.unit() - 0.5; norm += r * r; }
            norm = std::sqrt(norm);
            const double scale = 0.737998046875;  // evolve:param restart_scale real 0.01 1.0
            for (int i = 0; i < n; ++i) {
                double r = (rng.unit() - 0.5) / norm;
                h[i] = best[i] + scale * r;
                mass += h[i];
            }
            for (int i = 0; i < n; ++i) h[i] /= mass;
        }
        std::fill(velocity.begin(), velocity.end(), 0.0);
        // Draw a restart-specific rate from a log-uniform distribution.
        const double restart_rate = kRate * std::pow(kRateSpread, rng.unit() * 2.0 - 1.0);

        for (int step = 0; step < kSteps; ++step) {
            correlate(h, lags, c);
            const double worst = floor_of(c, lags);

            // Scored here rather than after the step, from the correlations
            // the gradient is about to use anyway. The vector every iteration
            // scores is the one the previous iteration produced, so nothing
            // goes unscored except the last, which is scored below the loop.
            {
                const double score = ratio(h, worst, in.per_unit);
                if (score > best_score) { best_score = score; best = h; }
            }

            // Softmin over the shifts, shifted by the minimum. The shift is
            // not a nicety: exp(-beta * c[k]) underflows at any useful beta,
            // and the ratios are what the gradient needs.
            const double beta =
                kBetaFrom * std::pow(kBetaTo / kBetaFrom,
                                     static_cast<double>(step) / static_cast<double>(kSteps));
            double total = 0.0;
            for (int k = 0; k <= lags; ++k) {
                soft[k] = std::exp(-beta * (c[k] - worst));
                total += soft[k];
            }
            for (int k = 0; k <= lags; ++k) soft[k] /= total;

            // d c[k] / d h[j] = h[j+k] + h[j-k], each term present only when
            // its index is inside the list, so the smoothed objective's
            // gradient is the softmin-weighted sum of those. ASCENT: the
            // objective is maximised here, which is the sign the sister
            // example has the other way.
            std::fill(gradient.begin(), gradient.end(), 0.0);
            for (int k = 0; k <= lags; ++k) {
                if (soft[k] < 1e-14) continue;
                const double weight = soft[k];
                for (int j = 0; j < n; ++j) {
                    if (j + k < n) gradient[j] += weight * h[j + k];
                    if (j - k >= 0) gradient[j] += weight * h[j - k];
                }
            }

            const double rate = restart_rate / (1.0 + 2e-5 * static_cast<double>(step));
            for (int i = 0; i < n; ++i) {
                velocity[i] = kMomentum * velocity[i] + rate * gradient[i];
                h[i] += velocity[i];
            }
            project(h, sorted);
        }

        // The vector the last step produced, which the loop never scored.
        correlate(h, lags, c);
        const double score = ratio(h, floor_of(c, lags), in.per_unit);
        if (score > best_score) { best_score = score; best = h; }
    }

    // OFF THE BOUNDARY OF CONDITION 3 BEFORE ANSWERING. The search works on the
    // unit simplex because that makes the step size mean something, but "sums
    // to one" and condition 3's "the sum is at least 1" are the same number,
    // and condition 3 has no tolerance — `problem.txt` says at least 1 and
    // means it. Measured on the sister example, where the same collision was
    // found: of four runs whose projection landed on 1, three came back
    // refused with sums of 0.9999999999999978, ...23 and ...97. The objective
    // is invariant under scale, so answering with twice the heights costs
    // exactly nothing and puts a whole factor between the answer and the line.
    double mass = 0.0;
    for (const double value : best) mass += value;
    const double scale = 2.0 / mass;
    for (double& value : best) value *= scale;

    // The claim is the floor this program computed, stated exactly, and taken
    // AFTER the rescale so that it describes the heights actually written.
    correlate(best, lags, c);
    const double claimed = floor_of(c, lags);

    std::ofstream output(argv[2]);
    output.precision(17);
    output << "{\"heights\":[";
    for (size_t i = 0; i < best.size(); ++i) {
        if (i) output << ',';
        output << best[i];
    }
    output << "],\"floor\":" << claimed << "}\n";
    return 0;
}
