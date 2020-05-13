% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% fourbarprisDE1
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% Jg_rot [3x1]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:10
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = fourbarprisDE1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),uint8(0),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE1_jacobig_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbarprisDE1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE1_jacobig_rot_sym_varpar: pkin has to be [3x1] (double)');
Jg_rot=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:19
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:20
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->12), mult. (11->10), div. (1->1), fcn. (2->2), ass. (0->4)
	t41 = (qJ(1) + pkin(3));
	t40 = (-pkin(2) - t41);
	t39 = (-pkin(2) + t41);
	t1 = [0; 0; (pkin(1) ^ 2 - pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3)) / t41 * (-(pkin(1) + t40) * (pkin(1) + t39) * (pkin(1) - t39) * (pkin(1) - t40)) ^ (-0.1e1 / 0.2e1);];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:20
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->12), mult. (11->10), div. (1->1), fcn. (2->2), ass. (0->4)
	t23 = (qJ(1) + pkin(3));
	t22 = (-pkin(2) - t23);
	t21 = (-pkin(2) + t23);
	t1 = [0; 0; (pkin(1) ^ 2 - pkin(2) ^ 2 - qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3)) / t23 * (-(pkin(1) + t22) * (pkin(1) + t21) * (pkin(1) - t21) * (pkin(1) - t22)) ^ (-0.1e1 / 0.2e1);];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 09:10:20
	% EndTime: 2020-05-07 09:10:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (14->14), mult. (6->6), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0; 0; -2 * (pkin(3) + qJ(1)) * (-(pkin(1) - pkin(2) - pkin(3) - qJ(1)) * (pkin(1) - pkin(2) + pkin(3) + qJ(1)) * (pkin(1) + pkin(2) - pkin(3) - qJ(1)) * (pkin(1) + pkin(2) + pkin(3) + qJ(1))) ^ (-0.1e1 / 0.2e1);];
	Jg_rot = t1;
end