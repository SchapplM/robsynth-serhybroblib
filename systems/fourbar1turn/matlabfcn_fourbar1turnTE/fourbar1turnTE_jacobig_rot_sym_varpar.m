% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% fourbar1turnTE
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% Jg_rot [3x2]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = fourbar1turnTE_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_jacobig_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnTE_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_jacobig_rot_sym_varpar: pkin has to be [5x1] (double)');
Jg_rot=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 1, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)); 0, -cos(qJ(1)); 1, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:26
	% EndTime: 2020-04-12 19:20:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (566->29), mult. (764->47), div. (28->6), fcn. (216->6), ass. (0->25)
	t92 = -2 * pkin(2);
	t91 = (-pkin(3) - pkin(4));
	t90 = (-pkin(3) + pkin(4));
	t80 = pkin(1) ^ 2;
	t76 = cos(qJ(2));
	t88 = pkin(1) * t76;
	t84 = t88 * t92 + t80;
	t71 = (pkin(2) ^ 2) + t84;
	t68 = pkin(3) ^ 2 - pkin(4) ^ 2 + t71;
	t72 = -pkin(2) + t88;
	t74 = sin(qJ(2));
	t66 = ((pkin(2) - t91) * (pkin(2) + t91)) + t84;
	t67 = ((pkin(2) - t90) * (pkin(2) + t90)) + t84;
	t81 = sqrt(-t66 * t67);
	t85 = t74 * t81;
	t60 = -pkin(1) * t85 - t72 * t68;
	t59 = 0.1e1 / t60 ^ 2;
	t87 = t74 * pkin(1);
	t65 = t68 * t87;
	t61 = -t72 * t81 + t65;
	t69 = 0.1e1 / t71;
	t83 = 0.1e1 / t71 ^ 2 * t92;
	t86 = 0.1e1 / t81 * (-t66 - t67) * pkin(2) * t87;
	t89 = (((0.2e1 * t80 * t74 ^ 2 * pkin(2) - t72 * t86) * t69 + ((t76 * t68 + t85) * t69 + t61 * t74 * t83) * pkin(1)) / t60 - (t65 * t69 + (-t81 * t76 * t69 + ((t72 * t92 - t86) * t69 + t60 * t83) * t74) * pkin(1)) * t61 * t59) / (t61 ^ 2 * t59 + 0.1e1) * t71 + 0.1e1;
	t1 = [0, t89 * sin(qJ(1)); 0, -t89 * cos(qJ(1)); 1, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:26
	% EndTime: 2020-04-12 19:20:26
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (565->28), mult. (768->47), div. (28->6), fcn. (214->6), ass. (0->25)
	t73 = -pkin(3) - pkin(4);
	t72 = -pkin(3) + pkin(4);
	t59 = sin(qJ(2));
	t71 = pkin(2) * t59;
	t60 = cos(qJ(2));
	t70 = pkin(2) * t60;
	t67 = (-0.2e1 * t70 + pkin(1)) * pkin(1);
	t51 = (pkin(2) - t73) * (pkin(2) + t73) + t67;
	t52 = (pkin(2) - t72) * (pkin(2) + t72) + t67;
	t64 = sqrt(-t51 * t52);
	t66 = pkin(1) * t71;
	t69 = (-t51 - t52) * t66 / t64;
	t68 = t59 * t64;
	t62 = pkin(2) ^ 2;
	t56 = t62 + t67;
	t57 = pkin(1) - t70;
	t55 = 0.1e1 / t56 ^ 2;
	t54 = 0.1e1 / t56;
	t53 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t56;
	t50 = t53 * t71;
	t46 = t57 * t64 + t50;
	t45 = -pkin(2) * t68 + t57 * t53;
	t44 = 0.1e1 / t45 ^ 2;
	t42 = 0.2e1 * (-((t57 * t69 + (t60 * t53 + t68) * pkin(2)) * t54 / 0.2e1 + (-pkin(2) * t46 * t55 + t54 * t62 * t59) * t59 * pkin(1)) / t45 - (-(t50 + (-t59 * t69 - t60 * t64) * pkin(2)) * t54 / 0.2e1 + (t45 * t55 - t54 * t57) * t66) * t46 * t44) / (t46 ^ 2 * t44 + 0.1e1) * t56;
	t1 = [0, t42 * sin(qJ(1)); 0, -t42 * cos(qJ(1)); 1, 0;];
	Jg_rot = t1;
end