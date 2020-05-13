% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh2m2OL
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh2m2OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh2m2OL_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2OL_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_jacobia_rot_sym_varpar: pkin has to be [5x1] (double)');
Ja_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:56
	% EndTime: 2020-05-03 06:34:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:58
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (1140->19), mult. (548->54), div. (133->9), fcn. (803->9), ass. (0->38)
	t65 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
	t64 = cos(t65);
	t63 = sin(t65);
	t67 = sin(qJ(1));
	t75 = t67 * t63;
	t59 = atan2(t75, t64);
	t52 = sin(t59);
	t53 = cos(t59);
	t49 = t52 * t75 + t53 * t64;
	t48 = 0.1e1 / t49 ^ 2;
	t69 = cos(qJ(1));
	t81 = t48 * t69 ^ 2;
	t80 = t52 * t64;
	t68 = cos(qJ(6));
	t71 = t69 * t68;
	t66 = sin(qJ(6));
	t74 = t67 * t66;
	t57 = t64 * t71 - t74;
	t55 = 0.1e1 / t57 ^ 2;
	t72 = t69 * t66;
	t73 = t67 * t68;
	t56 = t64 * t72 + t73;
	t79 = t55 * t56;
	t60 = t63 ^ 2;
	t78 = t60 / t64 ^ 2;
	t77 = t63 * t69;
	t58 = 0.1e1 / (t67 ^ 2 * t78 + 0.1e1);
	t76 = t67 * t58;
	t70 = t56 ^ 2 * t55 + 0.1e1;
	t61 = 0.1e1 / t64;
	t54 = 0.1e1 / t57;
	t51 = 0.1e1 / t70;
	t50 = (0.1e1 + t78) * t76;
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (t60 * t81 + 0.1e1);
	t45 = (-t54 * t66 + t68 * t79) * t51 * t77;
	t44 = (-t64 * t47 + (t67 * t80 - t53 * t63 + (t53 * t75 - t80) * t50) * t63 * t48) * t69 * t46;
	t1 = [t61 * t58 * t77, t50, t50, t50, t50, 0; (t47 * t75 + (t53 * t60 * t61 * t76 + (-t58 + 0.1e1) * t63 * t52) * t63 * t81) * t46, t44, t44, t44, t44, 0; ((-t64 * t74 + t71) * t54 - (-t64 * t73 - t72) * t79) * t51, t45, t45, t45, t45, t70 * t51;];
	Ja_rot = t1;
end