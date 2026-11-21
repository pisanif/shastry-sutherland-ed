%% ========================================================================
%  Exact diagonalisation (ED) of a finite Shastry-Sutherland cluster
%  -- fixed total-S^z sector (memory-efficient version) --
% =========================================================================
%
%  Reference : article about to be submitted
%
%  Companion of ED_SS_full_hilbert_space.m: same model and cluster, but the
%  Hamiltonian is built only in the sector of fixed total S^z.
%
%  MODEL
%  -----
%  Spin-1/2 Heisenberg model on a finite cluster of dimers (open boundaries)
%  of the Shastry-Sutherland lattice:
%
%      H =  J  * sum_{dimers (i,j)}  S_i . S_j     (intra-dimer bonds)
%         + J' * sum_{bonds  (i,j)}  S_i . S_j     (inter-dimer bonds)
%
%  The total S^z commutes with H, so H is block diagonal and only the block
%  with S^z = Sz_target is built (dimension C(N, N/2+Sz_target); for N = 24
%  and Sz_target = 1 this is 2,496,144 instead of 2^24 = 16,777,216). No
%  Zeeman term is needed, since the sector is selected explicitly.
%
%  WHAT THE SCRIPT DOES
%  --------------------
%   1. Sets parameters and the cluster geometry (bonds, and positions for N = 24).
%   2. Builds the basis of the S^z sector and the Hamiltonian (sparse).
%   3. Computes the lowest eigenstates with eigs.
%   4. (N = 24 only) Plots <S^z> of one selected state on the cluster.
%   5. Decomposes that state in the dimer basis (singlet/triplet on each
%      dimer), prints the dominant components and optionally saves them
%      to Excel files.
%   6. Computes the neutron-scattering intensity of the one-triplet part of
%      the state (ground state approximated by a product of dimer singlets).
%
%  CONVENTIONS
%  -----------
%   * Spins are numbered 1..N; spins (2m-1, 2m) form dimer m.
%   * A basis state is stored as an integer whose N binary digits encode the
%     spins: bit = 1 -> spin up, bit = 0 -> spin down. Spin 1 is the
%     most-significant bit, spin N the least-significant one.
%     (This is a different, bit-based indexing from the full-space script.)
%   * Dimer basis (per dimer): 0 = singlet s, 1 = triplet t+ (S^z=+1),
%     2 = triplet t0 (S^z=0), 3 = triplet t- (S^z=-1), with
%         s  = (|down,up> - |up,down>)/sqrt(2),   t0 = (|up,down> + |down,up>)/sqrt(2).
%     The sign of the singlet is opposite to the one used in the full-space
%     script, so amplitudes differ by (-1)^(number of singlets); the
%     intensities are not affected.
%   * Positions in Angstrom; the scattering vector is q = 2*pi/a * (h,k,0),
%     with h,k in reciprocal lattice units (r.l.u.).
%
%  NEUTRON CROSS SECTION (one-triplet components)
%  -------------------------------------------------
%  If the state is written as sum_j c_j |triplet on dimer j, singlets
%  elsewhere>, its overlap with sum_i exp(i q.r_i) S_i acting on a
%  product-of-singlets ground state gives
%       A(q) = sum_j c_j * 1/2 * ( exp(i q.r_{2j-1}) - exp(i q.r_{2j}) ),
%       I(q) = |A(q)|^2 .
%  Equivalently A(q) = i * sum_j c_j exp(i q.S_j) sin(q.d_j), with S_j the
%  centre of dimer j and d_j half the vector joining its two spins.
%  Polarisation factor, magnetic form factor and prefactors are not included.
%
%  REQUIREMENTS
%  ------------
%  MATLAB (base installation; no toolbox needed). Memory: a few GB for N = 24.
% =========================================================================

clear; close all; clc;

%% 1. PARAMETERS ----------------------------------------------------------
N         = 24;               % number of spins: 8 (small test), 24 or 32
Sz_target = 1;                % total S^z of the sector (integer, N even)

J       = 7.023;              % intra-dimer coupling
J_prime = 0.63*J;             % inter-dimer coupling (J'/J = 0.63)

num_eig      = 16;            % number of lowest eigenstates computed in the sector
target_state = 13;            % eigenstate analysed below (1 = lowest in the sector)

% Dimer-basis decomposition: keep components with |amplitude| > cutoff and
% at most max_triplets triplets.
cutoff_full       = 1e-4;     % general decomposition (listing/saving)
max_triplets_full = 4;
cutoff_1T         = 1e-5;     % one-triplet components used for the cross section

save_tables = false;          % true: write the decompositions to Excel files
output_dir  = pwd;            % folder for those files

a = 8.985;                    % in-plane lattice constant (Angstrom)

%% 2. CLUSTER GEOMETRY ----------------------------------------------------
n_dimers = N/2;

% Inter-dimer (J') bonds: each row is a pair of coupled spins [i, j].
switch N
    case 8
        J_prime_bonds = [ 1  3;  1  4;  3  5;  3  6;
                          5  7;  5  8;  1  7;  2  7];
    case 24
        J_prime_bonds = [10  1; 10  2;  8  1;  8  2;
                          2  3;  2  4; 12  3; 12  4;
                          4  5;  4  6; 14  5; 14  6;
                          6  8;  6  7; 16  8; 16  7;
                         17  9; 17 10;  3 17;  3 18;
                         19 12; 19 11;  5 19;  5 20;
                         21 14; 21 13;  7 21;  7 22;
                         23 16; 23 15;  1 23;  1 24];
    case 32
        J_prime_bonds = [10  1; 10  2;  7  1;  7  2;
                          2  3;  2  4; 11  3; 11  4;
                          4  5;  4  6; 13  5; 13  6;
                          6  7;  6  8; 15  7; 15  8;
                         17  9; 17 10;  3 17;  3 18;
                         19 11; 19 12;  5 19;  5 20;
                         21 13; 21 14;  8 21;  8 22;
                         23 15; 23 16;  1 23;  1 24;
                         18 25; 18 26; 26 11; 26 12;
                         20 27; 20 28; 28 13; 28 14;
                         22 29; 22 30; 29 15; 29 16;
                         24 31; 24 32; 32  9; 32 10];
    otherwise
        error('No bond table defined for N = %d.', N);
end

% Spin positions (x, y, z) in Angstrom, defined for the N = 24 cluster only
% (used for plotting and for the cross section). Rows 2m-1 and 2m are the
% two spins of dimer m.
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

%% 3. BASIS AND HAMILTONIAN -----------------------------------------------
fprintf('Building basis for Sz = %d ...\n', Sz_target);
basis = build_basis(N, Sz_target);        % sorted column vector of bit-encoded states
dim   = length(basis);
fprintf('Hilbert space dimension: %d\n', dim);

H = sparse(dim, dim);

% Intra-dimer terms: spins (1,2), (3,4), ...
for i = 1:2:N-1
    H = H + heisenberg_term(N, basis, i, i+1, J);
end

% Inter-dimer terms
for b = 1:size(J_prime_bonds, 1)
    H = H + heisenberg_term(N, basis, J_prime_bonds(b,1), J_prime_bonds(b,2), J_prime);
end

%% 4. DIAGONALISATION -----------------------------------------------------
% Lowest num_eig eigenstates of the sector (Lanczos). If convergence is a
% problem, eigs accepts a tolerance and a maximum number of iterations
% (e.g. opts.tol = 1e-11; opts.maxit = 1e4; eigs(H, num_eig, 'smallestreal', opts)).
fprintf('Diagonalizing ...\n');
tic
[V, D] = eigs(H, num_eig, 'smallestreal');
toc
[E, order] = sort(real(diag(D)));         % ascending energies
V = V(:, order);                          % V(:,n) is the n-th eigenstate of the sector

Egs=-3/4*N/2*J;   % ground state energy for the ss (exact product of singlets on the dimers)
fprintf('Lowest eigenvalues (state, E, E - Egs):\n');
fprintf(' %4d  %12.6f  %10.5f\n', [(1:num_eig).', E, E - Egs].');

if N ~= 24
    fprintf('Positions are only defined for N = 24: skipping the analysis.\n');
    return
end

state = V(:, target_state);               % state analysed below

%% 5. <S^z> OF THE SELECTED STATE -----------------------------------------
% S^z_i is diagonal in the basis: +1/2 for a 1-bit (up), -1/2 for a 0-bit.
% (<S^x_i> and <S^y_i> are identically zero in a state of definite S^z.)
prob = abs(state).^2;
exp_Sz = zeros(1, N);
for site = 1:N
    exp_Sz(site) = sum(prob .* (bitget(basis, N - site + 1) - 0.5));
end

% Plot: each spin is an arrow centred on its position, with length
% proportional to <S^z> (red: <S^z> > 0, blue: <S^z> <= 0). Dimers are
% labelled by their index; solid lines are dimer bonds, dashed lines are
% inter-dimer bonds.
arrow_scale = 10;                         % enlarges the arrows for readability

figure; hold on;
for k = 1:N
    dz = arrow_scale * exp_Sz(k);
    if dz > 0
        color = 'r';
    else
        color = 'b';
    end
    quiver(positions(k,1), positions(k,2) - dz/2, 0, dz, ...
           'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
end

for m = 1:n_dimers
    i = 2*m - 1;  j = 2*m;
    centre = 0.5 * (positions(i,1:2) + positions(j,1:2));
    text(centre(1) + 0.2, centre(2), num2str(m), 'FontSize', 15);   % dimer label
    plot(positions([i j], 1), positions([i j], 2), 'k-', 'LineWidth', 1.5);
end
for b = 1:size(J_prime_bonds, 1)
    i = J_prime_bonds(b,1);  j = J_prime_bonds(b,2);
    plot(positions([i j], 1), positions([i j], 2), 'k--');
end
axis equal;
axis off;
axis tight;
hold off;

%% 6. DIMER-BASIS DECOMPOSITION -------------------------------------------
% (a) general decomposition: dominant components with up to max_triplets triplets
fprintf('\nTransforming into the dimer basis (up to %d triplets)...\n', max_triplets_full);
tic
[configs, amps] = dimer_basis_decomposition(N, basis, state, cutoff_full, max_triplets_full);
toc
print_dimer_components(configs, amps, 30);
if save_tables
    save_dimer_components(configs, amps, fullfile(output_dir, ...
        sprintf('dimer_coeff_list_Sz_%d_%.0e_%dt_state_%d.xlsx', Sz_target, cutoff_full, max_triplets_full, target_state)));
end

% (b) one-triplet part of the state: one triplet on dimer j, singlets elsewhere
fprintf('\nTransforming into the dimer basis (one triplet only)...\n');
tic
[configs_1T, amps_1T] = dimer_basis_decomposition(N, basis, state, cutoff_1T, 1);
toc
print_dimer_components(configs_1T, amps_1T, 30);
if save_tables
    save_dimer_components(configs_1T, amps_1T, fullfile(output_dir, ...
        sprintf('dimer_coeff_list_Sz_%d_%.0e_1t_state_%d.xlsx', Sz_target, cutoff_1T, target_state)));
end

% c_1T(j): amplitude of the one-triplet component with the triplet on dimer j
% (components below cutoff_1T are taken as zero).
c_1T = zeros(1, n_dimers);
for n = 1:numel(configs_1T)
    j = find(configs_1T{n} ~= '0');       % dimer carrying the triplet
    if ~isempty(j)
        c_1T(j) = amps_1T(n);
    end
end

%% 7. NEUTRON CROSS SECTION, ONE-TRIPLET APPROXIMATION --------------------
h_vals = -4:0.05:4;                       % scattering-vector grid (r.l.u.)
k_vals = -4:0.05:4;
[h_grid, k_grid] = ndgrid(h_vals, k_vals);        % first index: h, second index: k
q = (2*pi/a) * [h_grid(:), k_grid(:), zeros(numel(h_grid), 1)];

r_first  = positions(1:2:end, :);         % first spin of each dimer  (n_dimers x 3)
r_second = positions(2:2:end, :);         % second spin of each dimer (n_dimers x 3)
amp_q = 0.5 * (exp(1i*q*r_first.') - exp(1i*q*r_second.')) * c_1T.';
CS = reshape(abs(amp_q).^2, numel(h_vals), numel(k_vals));      % CS(h index, k index)

% NOTE ON ORIENTATION: CS is indexed (h, k), but pcolor draws the first index
% along the vertical axis, so the figures below show h vertically and k
% horizontally (as in the original script). Use CS.' to swap them.
figure;
s = pcolor(h_vals, k_vals, CS);
s.EdgeColor = 'none';
colorbar;
axis equal tight;
title('Intensity, one-triplet approximation');

% Symmetrised pattern: CS + its 90-degree rotation, plus the left-right mirror
% of that sum (this is to account for the C4 symmetry of the lattice and
% the mirror-symmetry version of the pinwheel).
I_sym = CS + rot90(CS);
I_sym = I_sym + fliplr(I_sym);
figure;
s = pcolor(h_vals, k_vals, I_sym);
s.EdgeColor = 'none';
colorbar;
axis equal tight;
title('Intensity, one-triplet approximation (symmetrised)');

%% LOCAL FUNCTIONS --------------------------------------------------------

function basis = build_basis(N, Sz_target)
% All basis states of the sector S^z_tot = Sz_target, as a sorted column
% vector of integers (bit 1 = spin up, spin 1 = most-significant bit).
    n_up = N/2 + Sz_target;                    % number of up spins
    if mod(n_up, 1) ~= 0 || n_up < 0 || n_up > N
        error('Incompatible N and Sz_target');
    end
    combs = nchoosek(1:N, n_up);               % each row: positions of the up spins
    basis = sort(sum(2.^(N - combs), 2));      % set bit (N - position) for each up spin
end

function H = heisenberg_term(N, basis, i, j, J)
% Sparse matrix of J * S_i . S_j in the basis of the sector:
%   S_i.S_j = S^z_i S^z_j + 1/2 (S^+_i S^-_j + S^-_i S^+_j)
% The diagonal part is S^z_i S^z_j; the flip-flop part exchanges the spins
% of antiparallel pairs with matrix element J/2.
    dim   = length(basis);
    bit_i = bitget(basis, N - i + 1);          % 1 = up, 0 = down
    bit_j = bitget(basis, N - j + 1);

    diag_elements = J * (bit_i - 0.5) .* (bit_j - 0.5);

    rows    = find(bit_i ~= bit_j);            % states on which the flip-flop acts
    flipped = bitxor(bitxor(basis(rows), 2^(N-i)), 2^(N-j));
    [found, cols] = ismember(flipped, basis);  % column index of the flipped states
    rows = rows(found);
    cols = cols(found);

    H = sparse(1:dim, 1:dim, diag_elements, dim, dim) ...
      + sparse(rows, cols, 0.5*J, dim, dim);
end

function [configs, amps] = dimer_basis_decomposition(N, basis, evec, cutoff, max_triplets)
% Decomposes an eigenvector of the sector into the dimer basis.
%
%   evec          eigenvector in the spin basis (one column of V)
%   cutoff        spin-basis components with |amplitude| < cutoff are neglected
%   max_triplets  dimer states with more triplets than this are discarded
%
% Returns the dimer states as character strings configs{n} of N/2 digits
% (one per dimer, 0 = s, 1 = t+, 2 = t0, 3 = t-) and their complex
% amplitudes amps(n), sorted by decreasing |amplitude|.

    % Expansion of the spin state of a dimer ('s1 s2', 1 = up) in the dimer
    % basis: rows {dimer state, coefficient}
    spin_to_dimer = containers.Map;
    spin_to_dimer('11') = {[1, 1]};                           % t+
    spin_to_dimer('10') = {[2, 1/sqrt(2)], [0, -1/sqrt(2)]};  % t0 and s
    spin_to_dimer('01') = {[2, 1/sqrt(2)], [0,  1/sqrt(2)]};  % t0 and s
    spin_to_dimer('00') = {[3, 1]};                           % t-

    dimer_state_map = containers.Map;          % config string -> accumulated amplitude

    for idx = 1:length(evec)
        amp = evec(idx);
        if abs(amp) < cutoff
            continue;
        end
        spins = bitget(basis(idx), N:-1:1);    % spin string, spin 1 first (1 = up)

        % Expand the spin configuration dimer by dimer
        dimer_choices = {struct('config', [], 'amp', 1)};
        skip_state = false;
        for d = 1:2:N
            key = sprintf('%d%d', spins(d), spins(d+1));
            options = spin_to_dimer(key);

            new_choices = {};
            for k = 1:length(options)
                dimer_type = options{k}(1);
                coeff      = options{k}(2);
                for i = 1:length(dimer_choices)
                    item       = dimer_choices{i};
                    new_config = [item.config, dimer_type];
                    if sum(new_config ~= 0) > max_triplets
                        continue;
                    end
                    new_choices{end+1} = struct('config', new_config, 'amp', item.amp * coeff);
                end
            end
            dimer_choices = new_choices;

            if isempty(dimer_choices)
                skip_state = true;
                break;
            end
        end
        if skip_state
            continue;
        end

        % Accumulate the amplitude of every valid dimer configuration
        for i = 1:length(dimer_choices)
            key = sprintf('%d', dimer_choices{i}.config);
            amp_total = dimer_choices{i}.amp * amp;
            if isKey(dimer_state_map, key)
                dimer_state_map(key) = dimer_state_map(key) + amp_total;
            else
                dimer_state_map(key) = amp_total;
            end
        end
    end

    % Sort by decreasing |amplitude|
    configs = keys(dimer_state_map);
    amps    = cellfun(@(k) dimer_state_map(k), configs);
    [~, order] = sort(abs(amps), 'descend');
    configs = configs(order);
    amps    = amps(order);
end

function label = config_label(config)
% '010000000000' -> 's t+ s s s s s s s s s s'
    dimer_labels = {'s', 't+', 't0', 't-'};
    label = strjoin(dimer_labels((config - '0') + 1), ' ');
end

function print_dimer_components(configs, amps, n_max)
% Prints the (at most n_max) largest components of a dimer-basis decomposition.
    fprintf('\nMost relevant dimer basis states:\n');
    fprintf('  %-60s  %-25s  %-10s\n', 'Dimer State', 'Amplitude (complex)', '|Amplitude|^2');
    fprintf('  %s\n', repmat('-', 1, 100));
    for n = 1:min(n_max, numel(configs))
        c = amps(n);
        if abs(imag(c)) < 1e-10
            amp_str = sprintf('%+.5f', real(c));
        else
            amp_str = sprintf('%+.5f%+.5fi', real(c), imag(c));
        end
        fprintf('  %-60s  %-25s  %.5f\n', config_label(configs{n}), amp_str, abs(c)^2);
    end
end

function save_dimer_components(configs, amps, filename)
% Saves all components of a dimer-basis decomposition to a spreadsheet,
% sorted by decreasing |amplitude|.
    amps = amps(:);
    labels = cellfun(@config_label, configs(:), 'UniformOutput', false);
    T = table(labels, real(amps), imag(amps), abs(amps).^2, ...
        'VariableNames', {'DimerState', 'ReAmplitude', 'ImAmplitude', 'AbsAmplitudeSq'});
    writetable(T, filename);
    fprintf('\nResults saved to %s\n', filename);
end
