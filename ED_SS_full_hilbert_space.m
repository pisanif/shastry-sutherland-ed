%% ========================================================================
%  Exact diagonalisation (ED) of a finite Shastry-Sutherland cluster
%  -- full Hilbert space (2^24 states), no symmetry used --
% =========================================================================
%
%  Reference : Article about to be submitted
%
%  MODEL
%  -----
%  Spin-1/2 Heisenberg model on a finite cluster of 12 orthogonal dimers
%  (24 spins) of the Shastry-Sutherland lattice (open boundaries):
%
%      H =  J  * sum_{dimers (i,j)}  S_i . S_j     (intra-dimer bonds)
%         + J' * sum_{bonds  (i,j)}  S_i . S_j     (inter-dimer bonds)
%         -      sum_i  B . S_i                    (Zeeman term, g*mu_B absorbed in B)
%
%  The geometry (a = 8.985 Angstrom, Cu-Cu dimer distance 2.9 Angstrom)
%  corresponds to SrCu2(BO3)2.
%
%  WHAT THE SCRIPT DOES
%  --------------------
%   1. Defines parameters and the cluster geometry (positions, bonds).
%   2. Builds the single-spin operators and the Hamiltonian as sparse
%      matrices in the full 2^24-dimensional Hilbert space.
%   3. Computes the lowest eigenstates with eigs and prints their energies
%      and total S^z.
%   4. Plots the spin expectation values <S^x>, <S^z> of one selected state.
%   5. Rotates that state into the dimer basis (singlet/triplet on each
%      dimer) and lists its dominant components.
%   6. Computes the magnetic neutron-scattering intensity of the transition
%      ground state -> selected state:
%        (a) fast way: only the one-triplet components of the selected
%            state, with the ground state taken as a product of dimer singlets;
%        (b) exact: full ED wavefunctions of both states.
%
%  CONVENTIONS
%  -----------
%   * Spins are numbered 1..24; spins (2m-1, 2m) form dimer m.
%   * Single-spin basis: |up> = (1,0), |down> = (0,1)  (S^z = +1/2, -1/2).
%     Spin 1 is the left-most Kronecker factor: the basis-state index n of
%     the full space is read from the binary digits of (n-1), left-most
%     digit = spin 1, digit 0 = up, digit 1 = down.
%   * Dimer basis (per dimer): 0 = singlet, 1 = triplet t+ (S^z=+1),
%     2 = triplet t0 (S^z=0), 3 = triplet t- (S^z=-1). The index of the
%     12-dimer basis state is read from the base-4 digits of (n-1), with
%     dimer 1 as left-most digit.
%   * Energies are in the units of J; positions in Angstrom; the scattering
%     vector q = 2*pi/a * (h,k,0) is in-plane, with h,k in reciprocal
%     lattice units (r.l.u.).
%
%  NEUTRON CROSS SECTION
%  ---------------------
%  The intensity plotted is
%       I(q) = sum_{alpha=x,y,z} | < f | sum_i exp(i q.r_i) S_i^alpha | 0 > |^2 ,
%  with |0> the ground state and |f> the selected excited state. The
%  polarisation factor, magnetic form factor and energy-dependent prefactors
%  are NOT included.
%
%  REQUIREMENTS
%  ------------
%  MATLAB (base installation; no toolbox needed). Memory is the limiting
%  factor: the 72 sparse single-spin operators alone need several tens of GB
%  of RAM. The memory-efficient version of this calculation, restricted to a
%  fixed-S^z sector, is provided in the companion script.
% =========================================================================

clear; close all; clc;

%% 1. PARAMETERS ----------------------------------------------------------
n_spins  = 24;                % number of spins; Hilbert-space dimension is 2^n_spins
n_dimers = n_spins/2;         % number of dimers

J       = 7.023;              % intra-dimer coupling
J_prime = 0.63*J;             % inter-dimer coupling (J'/J = 0.63)

% Magnetic field B = [Bx, By, Bz]. A tiny field along z lifts the (2S+1)-fold
% degeneracy of the multiplets, so that every computed eigenstate has a
% definite total S^z (with H = ... - B.S and Bz > 0, S^z = +1 is lowest).
B = [0, 0, 1e-11];

n_states     = 16;            % number of lowest eigenstates to compute
target_state = 2;             % eigenstate to analyse (1 = ground state, 2 = next, ...)

a = 8.985;                    % in-plane lattice constant (Angstrom)

%% 2. CLUSTER GEOMETRY ----------------------------------------------------
% Spin positions (x, y, z) in Angstrom for SrCu2(BO3)2. Rows 2m-1 and 2m are the two spins
% of dimer m. The z component is set to zero because it is not relevant for the
% calcuclation  here.
positions = [-5.5179,   -1.0254,    0;    % dimer 1
             -3.4671,    1.0254,    0;
             -1.0254,    5.5179,    0;    % dimer 2
              1.0254,    3.4671,    0;
              5.51787,   1.02537,   0;    % dimer 3
              3.46713,  -1.02537,   0;
              1.02537,  -5.51787,   0;    % dimer 4
             -1.02537,  -3.46713,   0;
            -10.0104,    5.51787,   0;    % dimer 5
             -7.95963,   3.46713,   0;
              5.51787,  10.0104,    0;    % dimer 6
              3.46713,   7.95963,   0;
             10.0104,   -5.51787,   0;    % dimer 7
              7.95963,  -3.46713,   0;
             -5.51787, -10.0104,    0;    % dimer 8
             -3.46713,  -7.95963,   0;
             -5.51787,   7.95963,   0;    % dimer 9
             -3.46713,  10.0104,    0;
              7.95963,   5.51787,   0;    % dimer 10
             10.0104,    3.46713,   0;
              5.51787,  -7.95963,   0;    % dimer 11
              3.46713, -10.0104,    0;
             -7.95963,  -5.51787,   0;    % dimer 12
            -10.0104,   -3.46713,   0];

% Inter-dimer (J') bonds: each row is a pair of coupled spins [i, j].
J_prime_bonds = [10  1;  10  2;   8  1;   8  2;
                  2  3;   2  4;  12  3;  12  4;
                  4  5;   4  6;  14  5;  14  6;
                  6  8;   6  7;  16  8;  16  7;
                 17  9;  17 10;   3 17;   3 18;
                 19 12;  19 11;   5 19;   5 20;
                 21 14;  21 13;   7 21;   7 22;
                 23 16;  23 15;   1 23;   1 24];

%% 3. SPIN OPERATORS AND HAMILTONIAN --------------------------------------
% Single-spin-1/2 operators
Sx = 0.5*sparse([0 1; 1 0]);
Sy = 0.5*sparse([0 -1i; 1i 0]);
Sz = 0.5*sparse([1 0; 0 -1]);

% S_ops{k,1:3} = S^x, S^y, S^z of spin k, embedded in the full Hilbert space
fprintf('Building spin operators...\n');
S_ops = cell(n_spins, 3);
for k = 1:n_spins
    S_ops{k,1} = spin_operator(k, Sx, n_spins);
    S_ops{k,2} = spin_operator(k, Sy, n_spins);
    S_ops{k,3} = spin_operator(k, Sz, n_spins);
end

% Scalar product S_i . S_j of two spins
S_dot = @(i, j) S_ops{i,1}*S_ops{j,1} + S_ops{i,2}*S_ops{j,2} + S_ops{i,3}*S_ops{j,3};

fprintf('Building Hamiltonian...\n');
H = sparse(2^n_spins, 2^n_spins);

% Intra-dimer terms: spins (1,2), (3,4), ...
for m = 1:n_dimers
    H = H + J * S_dot(2*m-1, 2*m);
end

% Zeeman term (only the non-zero components of B)
for k = 1:n_spins
    for c = find(B ~= 0)
        H = H - B(c) * S_ops{k,c};
    end
end

% Inter-dimer terms
for b = 1:size(J_prime_bonds, 1)
    H = H + J_prime * S_dot(J_prime_bonds(b,1), J_prime_bonds(b,2));
end

%% 4. DIAGONALISATION -----------------------------------------------------
% Lowest n_states eigenstates (Lanczos). 'smallestreal' makes sure that the
% lowest-energy states are returned (eigs' default would be the largest
% magnitude). Increase n_states if more levels are needed.
fprintf('Diagonalising...\n');
tic
[V, D] = eigs(H, n_states, 'smallestreal');
toc
[E, order] = sort(real(diag(D)));    % ascending energies
V = V(:, order);                     % V(:,n) is the n-th eigenstate

% Overview of the spectrum: energy relative to the ground state and total S^z
Sz_total_diag = zeros(2^n_spins, 1);  % S^z_tot is diagonal in this basis
for k = 1:n_spins
    Sz_total_diag = Sz_total_diag + full(diag(S_ops{k,3}));
end

fprintf('\n state    E - E0      <Sz_tot>\n');
for n = 1:n_states
    Sz_n = real(sum(abs(V(:,n)).^2 .* Sz_total_diag));
    fprintf(' %4d  %10.5f  %9.4f\n', n, E(n) - E(1), Sz_n);
end

ground_state = V(:, 1);
state        = V(:, target_state);   % excited state analysed below

%% 5. SPIN EXPECTATION VALUES OF THE SELECTED STATE -----------------------
exp_Sx = zeros(1, n_spins);
exp_Sz = zeros(1, n_spins);
for k = 1:n_spins
    exp_Sx(k) = real(state' * (S_ops{k,1} * state));
    exp_Sz(k) = real(state' * (S_ops{k,3} * state));
end

% Plot: each spin is drawn as an arrow centred on its position, pointing
% along (<S^x>, <S^z>) (S^x along the horizontal, S^z along the vertical
% plot axis). Red: <S^z> > 0, blue: <S^z> <= 0. Solid lines are dimer bonds,
% dashed lines are inter-dimer bonds.
arrow_scale = 20;                    % enlarges the arrows for readability

figure; hold on;
for k = 1:n_spins
    dx = arrow_scale * exp_Sx(k);
    dz = arrow_scale * exp_Sz(k);
    if dz > 0
        color = 'r';
    else
        color = 'b';
    end
    quiver(positions(k,1) - dx/2, positions(k,2) - dz/2, dx, dz, ...
           'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
    text(positions(k,1) + 0.4, positions(k,2), num2str(k));   % spin label
end

for m = 1:n_dimers                   % dimer bonds (solid)
    i = 2*m - 1;  j = 2*m;
    plot(positions([i j], 1), positions([i j], 2), 'k-', 'LineWidth', 1.5);
end
for b = 1:size(J_prime_bonds, 1)     % inter-dimer bonds (dashed)
    i = J_prime_bonds(b,1);  j = J_prime_bonds(b,2);
    plot(positions([i j], 1), positions([i j], 2), 'k--');
end
xlabel('x [Angstrom]');
ylabel('y [Angstrom]');
title(sprintf('Spin expectation values, state %d', target_state));
axis equal;
hold off;

%% 6. DIMER BASIS ---------------------------------------------------------
% Local change of basis on every dimer, from {|uu>, |ud>, |du>, |dd>} to
% {singlet, t+, t0, t-}. Row r of T is the r-th new basis vector, i.e.
% new_coefficient(r) = sum_c T(r,c) * old_coefficient(c).
T = [0   1/sqrt(2)  -1/sqrt(2)   0;     % singlet (|ud> - |du>)/sqrt(2)
     1   0            0          0;     % t+      |uu>
     0   1/sqrt(2)   1/sqrt(2)   0;     % t0      (|ud> + |du>)/sqrt(2)
     0   0            0          1];    % t-      |dd>

% Same as multiplying by kron(T, T, ..., T) (12 factors), but without ever
% building that huge matrix.
state_dm = to_dimer_basis(state, T, n_dimers);

% List the dominant components in the dimer basis
threshold = 0.1;                     % only components with |amplitude| > threshold
labels    = 's+0-';                  % character used for each local dimer state
idx = find(abs(state_dm) > threshold);
[~, o] = sort(abs(state_dm(idx)), 'descend');
idx = idx(o);

fprintf('\nDominant components of state %d in the dimer basis (dimers 1..12):\n', target_state);
fprintf('   amplitude    weight   dimer states (s=singlet, +/0/-=triplet)\n');
for n = 1:numel(idx)
    digits = dec2base(idx(n) - 1, 4, n_dimers) - '0';    % local state of each dimer
    fprintf('  %9.4f  %8.4f   %s\n', real(state_dm(idx(n))), abs(state_dm(idx(n)))^2, labels(digits + 1));
end

%% 7. NEUTRON CROSS SECTION, ONE-TRIPLET APPROXIMATION --------------------
% One-triplet basis state: dimer m in triplet t (t = 1,2,3 for t+, t0, t-),
% all other dimers in the singlet. Its index in the dimer basis is
%       t * 4^(n_dimers - m) + 1.
amp_1T = zeros(3, n_dimers);         % amp_1T(t,m): amplitude of that basis state
for t = 1:3
    for m = 1:n_dimers
        amp_1T(t,m) = state_dm(t * 4^(n_dimers - m) + 1);
    end
end
fprintf('\nTotal weight of one-triplet components: %.4f\n', sum(abs(amp_1T(:)).^2));

% The selected state is a triplet with a definite S^z, so only one flavour
% (t+, t0 or t-) has non-zero amplitudes: pick the one with the largest weight.
[~, t_main] = max(sum(abs(amp_1T).^2, 2));
c_1T = amp_1T(t_main, :);            % one-triplet amplitude on each dimer (1 x 12)

% Matrix element of sum_i exp(i q.r_i) S_i between a singlet-product ground
% state and a triplet localised on dimer m: (exp(i q.r_1) - exp(i q.r_2))/2,
% with r_1, r_2 the positions of the two spins of the dimer. The amplitudes
% of the dimers are added coherently with the weights c_1T.
h_vals = -4:0.1:4;                   % scattering-vector grid (r.l.u.)
k_vals = -4:0.1:4;
[h_grid, k_grid] = meshgrid(h_vals, k_vals);
q = (2*pi/a) * [h_grid(:), k_grid(:), zeros(numel(h_grid), 1)];

r_first  = positions(1:2:end, :);    % first spin of each dimer  (12 x 3)
r_second = positions(2:2:end, :);    % second spin of each dimer (12 x 3)
amp_q = 0.5 * (exp(1i*q*r_first.') - exp(1i*q*r_second.')) * c_1T.';
CS_1T = reshape(abs(amp_q).^2, numel(k_vals), numel(h_vals));   % rows: k, columns: h

CS_1T_tot = CS_1T + rot90(CS_1T) + flipud(CS_1T + rot90(CS_1T));
figure;
s = pcolor(h_vals, k_vals, CS_1T_tot);
s.EdgeColor = 'none';
colorbar;
axis equal tight;
xlabel('h (r.l.u.)');
ylabel('k (r.l.u.)');
title('Intensity, one-triplet approximation');

%% 8. NEUTRON CROSS SECTION, FULL ED WAVEFUNCTIONS ------------------------
% I(q) = sum_alpha | sum_i exp(i q.r_i) <state| S_i^alpha |ground_state> |^2.
% The matrix elements <state|S_i^alpha|ground_state> do not depend on q, so
% they are computed once and the q-dependence is just a sum of phase factors.
fprintf('\nComputing matrix elements for the scattering cross section...\n');
M = zeros(n_spins, 3);
for k = 1:n_spins
    for alpha = 1:3
        M(k, alpha) = state' * (S_ops{k,alpha} * ground_state);
    end
end

h_vals = -4:0.25:4;                  % coarser grid (r.l.u.)
k_vals = -4:0.25:4;
[h_grid, k_grid] = meshgrid(h_vals, k_vals);
q = (2*pi/a) * [h_grid(:), k_grid(:), zeros(numel(h_grid), 1)];

phase_factors = exp(1i * (q * positions.'));       % (n_q x n_spins): exp(i q.r_i)
IN = sum(abs(phase_factors * M).^2, 2);            % sum over alpha = x, y, z
IN = reshape(IN, numel(k_vals), numel(h_vals));    % rows: k, columns: h

% Impose tzhe C4 lattice symmetry and symmetrise with the mirror-symmetric pattern (it accounts for the two
% possible pinwheels in the real lattice)

IN_sym = IN + rot90(IN) + flipud(IN + rot90(IN));

figure;
s = pcolor(h_vals, k_vals, IN_sym);
s.EdgeColor = 'none';
colorbar;
axis equal tight;
xlabel('h (r.l.u.)');
ylabel('k (r.l.u.)');
title('Intensity, full ED wavefunctions');

%% LOCAL FUNCTIONS --------------------------------------------------------

function op_full = spin_operator(spin_index, op, n_spins)
% Embeds the single-spin operator OP (2x2) acting on spin SPIN_INDEX into
% the full Hilbert space: 1 (x) ... (x) OP (x) ... (x) 1 (Kronecker products).
    op_full = kron(kron(speye(2^(spin_index-1)), op), speye(2^(n_spins-spin_index)));
end

function psi = to_dimer_basis(psi, T, n_dimers)
% Applies the 4x4 matrix T to every dimer of the state vector PSI (length
% 4^n_dimers), i.e. computes kron(T,...,T)*psi without building the matrix.
% Each dimer corresponds to one base-4 digit of the basis-state index.
    for m = 1:n_dimers
        psi = reshape(psi, 4^(m-1), 4, []);    % (lower digits, digit m, upper digits)
        psi = permute(psi, [2 1 3]);           % (digit m, lower digits, upper digits)
        psi = T * reshape(psi, 4, []);         % apply T on digit m
        psi = reshape(psi, 4, 4^(m-1), []);
        psi = permute(psi, [2 1 3]);           % back to the original ordering
    end
    psi = psi(:);
end
