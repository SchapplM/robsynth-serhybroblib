% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% fourbar1turnDE2
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = fourbar1turnDE2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_jacobig_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnDE2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_jacobig_rot_sym_varpar: pkin has to be [5x1] (double)');
Jg_rot=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 1, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(qJ(1)); 0, -cos(qJ(1)); 1, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:43
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (566->29), mult. (764->47), div. (28->6), fcn. (216->6), ass. (0->25)
	t78 = -2 * pkin(2);
	t77 = (-pkin(3) - pkin(4));
	t76 = (-pkin(3) + pkin(4));
	t66 = pkin(1) ^ 2;
	t62 = cos(qJ(2));
	t73 = t62 * pkin(1);
	t70 = t73 * t78 + t66;
	t57 = (pkin(2) ^ 2) + t70;
	t54 = pkin(3) ^ 2 - pkin(4) ^ 2 + t57;
	t58 = -pkin(2) + t73;
	t60 = sin(qJ(2));
	t52 = ((pkin(2) - t77) * (pkin(2) + t77)) + t70;
	t53 = ((pkin(2) - t76) * (pkin(2) + t76)) + t70;
	t67 = sqrt(-t52 * t53);
	t71 = t60 * t67;
	t46 = -pkin(1) * t71 - t54 * t58;
	t45 = 0.1e1 / t46 ^ 2;
	t74 = t60 * pkin(1);
	t51 = t54 * t74;
	t47 = -t58 * t67 + t51;
	t55 = 0.1e1 / t57;
	t69 = 0.1e1 / t57 ^ 2 * t78;
	t72 = 0.1e1 / t67 * (-t52 - t53) * pkin(2) * t74;
	t75 = (((0.2e1 * t66 * t60 ^ 2 * pkin(2) - t58 * t72) * t55 + ((t62 * t54 + t71) * t55 + t47 * t60 * t69) * pkin(1)) / t46 - (t51 * t55 + (-t67 * t62 * t55 + ((t58 * t78 - t72) * t55 + t46 * t69) * t60) * pkin(1)) * t47 * t45) / (t45 * t47 ^ 2 + 0.1e1) * t57 + 0.1e1;
	t1 = [0, t75 * sin(qJ(1)); 0, -t75 * cos(qJ(1)); 1, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:43
	% EndTime: 2020-04-12 19:35:43
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (565->28), mult. (768->47), div. (28->6), fcn. (214->6), ass. (0->25)
	t89 = -pkin(3) - pkin(4);
	t88 = -pkin(3) + pkin(4);
	t75 = sin(qJ(2));
	t87 = pkin(2) * t75;
	t76 = cos(qJ(2));
	t86 = pkin(2) * t76;
	t83 = (-0.2e1 * t86 + pkin(1)) * pkin(1);
	t67 = (pkin(2) - t89) * (pkin(2) + t89) + t83;
	t68 = (pkin(2) - t88) * (pkin(2) + t88) + t83;
	t80 = sqrt(-t67 * t68);
	t82 = pkin(1) * t87;
	t85 = (-t67 - t68) * t82 / t80;
	t84 = t75 * t80;
	t78 = pkin(2) ^ 2;
	t72 = t78 + t83;
	t73 = pkin(1) - t86;
	t71 = 0.1e1 / t72 ^ 2;
	t70 = 0.1e1 / t72;
	t69 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t72;
	t66 = t69 * t87;
	t62 = t73 * t80 + t66;
	t61 = -pkin(2) * t84 + t73 * t69;
	t60 = 0.1e1 / t61 ^ 2;
	t58 = 0.2e1 * (-((t73 * t85 + (t76 * t69 + t84) * pkin(2)) * t70 / 0.2e1 + (-pkin(2) * t62 * t71 + t70 * t78 * t75) * t75 * pkin(1)) / t61 - (-(t66 + (-t75 * t85 - t76 * t80) * pkin(2)) * t70 / 0.2e1 + (t61 * t71 - t70 * t73) * t82) * t62 * t60) / (t62 ^ 2 * t60 + 0.1e1) * t72;
	t1 = [0, t58 * sin(qJ(1)); 0, -t58 * cos(qJ(1)); 1, 0;];
	Jg_rot = t1;
end