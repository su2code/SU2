#!/usr/bin/env python3
"""Step 4: train the network offline on the inverted beta and write nn_weights.dat.

One SU2 iteration from restart_inverted.csv writes the four features at every point (features.csv). The points
that FILTER_SHIELD keeps are split into training and validation sets. The scaling (Box-Cox lambda by maximum
likelihood, then z-score), the network (tanh layers, beta = 1 + tanh(output)) and the file layout are those of
CTurbSASolver, so SU2 reproduces the predicted beta.
"""

import os

import numpy as np

from fiml_case import read_restart, require, run_parser, run_su2

FEATURES = ["FIML_Prod_Dest", "FIML_Chi", "FIML_Delta", "FIML_Strain_Vort"]
BOX_COX_FLOOR = 1e-12


def box_cox(x, lam):
    x = np.maximum(x, BOX_COX_FLOOR)
    return np.log(x) if abs(lam) < 1e-10 else (x**lam - 1.0) / lam


def fit_scaling(x):
    """Box-Cox lambda on the grid used by SU2, then the mean and std of the transformed feature."""
    log_sum = np.log(np.maximum(x, BOX_COX_FLOOR)).sum()
    best = None
    for lam in np.linspace(-2.0, 2.0, 81):
        y = box_cox(x, lam)
        var = max(y.var(), 1e-300)
        ll = -0.5 * x.size * np.log(var) + (lam - 1.0) * log_sum
        if best is None or ll > best[0]:
            best = (ll, lam, y.mean(), max(np.sqrt(var), 1e-12))
    return best[1:]


class Network:
    """Fully connected tanh network; W[l] has the bias as its last row, as in nn_weights.dat."""

    def __init__(self, n_in, n_neurons, n_hidden, rng):
        sizes = [n_in] + [n_neurons] * n_hidden
        self.W = []
        for a, b in zip(sizes[:-1], sizes[1:]):
            lim = np.sqrt(6.0 / (a + b))
            self.W.append(np.vstack([rng.uniform(-lim, lim, (a, b)), np.zeros((1, b))]))
        self.W.append(np.zeros((n_neurons + 1, 1)))

    def forward(self, x):
        acts = [x]
        for k, W in enumerate(self.W):
            z = np.hstack([acts[-1], np.ones((len(x), 1))]) @ W
            acts.append(np.tanh(z))
        return 1.0 + acts[-1][:, 0], acts

    def gradient(self, x, target):
        beta, acts = self.forward(x)
        delta = ((beta - target) * (1.0 - acts[-1][:, 0] ** 2))[:, None] / len(x)
        grads = [None] * len(self.W)
        for k in range(len(self.W) - 1, -1, -1):
            a = np.hstack([acts[k], np.ones((len(x), 1))])
            grads[k] = a.T @ delta
            if k > 0:
                delta = (delta @ self.W[k][:-1].T) * (1.0 - acts[k] ** 2)
        return 0.5 * np.mean((beta - target) ** 2), grads


def train(net, x, y, x_val, y_val, epochs, batch, lr, rng):
    m = [np.zeros_like(W) for W in net.W]
    v = [np.zeros_like(W) for W in net.W]
    step = 0
    for epoch in range(epochs):
        order = rng.permutation(len(x))
        for start in range(0, len(x), batch):
            idx = order[start:start + batch]
            _, grads = net.gradient(x[idx], y[idx])
            step += 1
            for k, g in enumerate(grads):
                m[k] = 0.9 * m[k] + 0.1 * g
                v[k] = 0.999 * v[k] + 0.001 * g * g
                net.W[k] -= lr * (m[k] / (1 - 0.9**step)) / (np.sqrt(v[k] / (1 - 0.999**step)) + 1e-8)
        if epoch % 50 == 0 or epoch == epochs - 1:
            mse = np.mean((net.forward(x)[0] - y) ** 2)
            mse_val = np.mean((net.forward(x_val)[0] - y_val) ** 2)
            print(f"epoch {epoch:4d}: training MSE {mse:.4e}, validation MSE {mse_val:.4e}", flush=True)


def write_weights(path, net, mean, std, lam, n_neurons):
    with open(path, "w") as f:
        f.write(f"{len(mean)} Inputs,{len(net.W) - 1} Hidden Layers,{n_neurons} Hidden Nodes\n")
        f.write(" ".join(f"{v:.15g}" for v in [*mean, *std, *lam]) + "\n")
        for W in net.W:
            f.writelines(f"{w:.15e}\n" for w in W.ravel())


def main():
    p = run_parser(__doc__)
    p.add_argument("--neurons", type=int, default=20)
    p.add_argument("--hidden", type=int, default=3)
    p.add_argument("--epochs", type=int, default=400)
    p.add_argument("--batch", type=int, default=512)
    p.add_argument("--lr", type=float, default=2e-3)
    p.add_argument("--validation", type=float, default=0.2, help="fraction of points held out")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    if not os.path.exists("features.csv"):
        require("restart_inverted.csv")
        run_su2(args, "features", {"SOLUTION_FILENAME": "restart_inverted", "RESTART_FILENAME": "features",
                                   "ITER": "1", "WRT_RESTART_COMPACT": "NO",
                                   "VOLUME_OUTPUT": "COORDINATES, SOLUTION, FIML, FIML_FEATURES",
                                   "OUTPUT_FILES": "RESTART_ASCII"})
    data = read_restart("features.csv")
    feats = np.column_stack([data[k] for k in FEATURES])
    beta = data["Beta_FIML"]

    shield = (feats[:, 2] >= 0.25) & (feats[:, 3] <= 2.0)
    outside = ~shield & (abs(beta - 1.0) > 0.05)
    print(f"{shield.sum()} of {beta.size} points inside the filter; "
          f"{outside.sum()} points outside it have |beta - 1| > 0.05 and stay at beta = 1")

    scaling = [fit_scaling(feats[shield, k]) for k in range(4)]
    lam, mean, std = (np.array(v) for v in zip(*scaling))
    x = np.column_stack([(box_cox(feats[shield, k], lam[k]) - mean[k]) / std[k] for k in range(4)])
    y = beta[shield]
    print("Box-Cox lambda", np.round(lam, 3), "mean", np.round(mean, 3), "std", np.round(std, 3))

    rng = np.random.default_rng(args.seed)
    val = rng.random(len(y)) < args.validation
    net = Network(4, args.neurons, args.hidden, rng)
    train(net, x[~val], y[~val], x[val], y[val], args.epochs, args.batch, args.lr, rng)
    write_weights("nn_weights.dat", net, mean, std, lam, args.neurons)

    pred = net.forward(x)[0]
    r2 = 1.0 - np.sum((pred - y) ** 2) / np.sum((y - y.mean()) ** 2)
    print(f"inside the filter: RMS(beta_nn - beta_inverted) = {np.sqrt(np.mean((pred - y) ** 2)):.4f}, "
          f"R^2 = {r2:.3f}, max error = {abs(pred - y).max():.3f}")
    print("wrote nn_weights.dat")


if __name__ == "__main__":
    main()
