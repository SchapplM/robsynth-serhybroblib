% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m2TE
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh3m2TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh3m2TE_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2TE_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_jacobia_rot_sym_varpar: pkin has to be [18x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:40
	% EndTime: 2020-05-07 01:41:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:40
	% EndTime: 2020-05-07 01:41:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:42
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:43
	% EndTime: 2020-05-07 01:41:43
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (356->19), mult. (520->43), div. (25->9), fcn. (765->13), ass. (0->33)
	t63 = sin(pkin(16));
	t64 = cos(pkin(16));
	t67 = sin(pkin(15));
	t70 = cos(pkin(15));
	t58 = t63 * t70 + t64 * t67;
	t59 = -t63 * t67 + t64 * t70;
	t62 = pkin(17) + pkin(18);
	t60 = sin(t62);
	t61 = cos(t62);
	t81 = -t58 * t60 + t61 * t59;
	t56 = t58 * t61 + t59 * t60;
	t66 = sin(qJ(1));
	t76 = t66 * t56;
	t47 = atan2(-t76, t81);
	t45 = sin(t47);
	t46 = cos(t47);
	t43 = -t45 * t76 + t46 * t81;
	t69 = cos(qJ(1));
	t80 = 0.1e1 / t43 ^ 2 * t69 ^ 2;
	t54 = t56 ^ 2;
	t48 = 0.1e1 / (0.1e1 + t66 ^ 2 * t54 / t81 ^ 2);
	t79 = t48 / t81;
	t65 = sin(qJ(4));
	t75 = t69 * t65;
	t68 = cos(qJ(4));
	t74 = t69 * t68;
	t52 = t66 * t65 - t74 * t81;
	t49 = 0.1e1 / t52 ^ 2;
	t50 = t66 * t68 + t81 * t75;
	t73 = t50 ^ 2 * t49 + 0.1e1;
	t71 = t81 * t66;
	t44 = 0.1e1 / t73;
	t1 = [-t69 * t56 * t79, 0, 0, 0; (-0.1e1 / t43 * t76 - (t46 * t54 * t66 * t79 + (t48 - 0.1e1) * t56 * t45) * t56 * t80) / (t54 * t80 + 0.1e1), 0, 0, 0; ((t65 * t71 - t74) / t52 + (t68 * t71 + t75) * t50 * t49) * t44, 0, 0, t73 * t44;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:44
	% EndTime: 2020-05-07 01:41:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end