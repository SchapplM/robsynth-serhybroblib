% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh2m1OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh2m1OL_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1OL_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_jacobia_rot_sym_varpar: pkin has to be [6x1] (double)');
Ja_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (695->21), mult. (440->54), div. (106->9), fcn. (646->9), ass. (0->38)
	t72 = qJ(2) + qJ(3) + qJ(4);
	t71 = cos(t72);
	t70 = sin(t72);
	t74 = sin(qJ(1));
	t82 = t74 * t70;
	t65 = atan2(-t82, t71);
	t61 = sin(t65);
	t62 = cos(t65);
	t56 = -t61 * t82 + t62 * t71;
	t55 = 0.1e1 / t56 ^ 2;
	t76 = cos(qJ(1));
	t88 = t55 * t76 ^ 2;
	t75 = cos(qJ(5));
	t78 = t76 * t75;
	t73 = sin(qJ(5));
	t81 = t74 * t73;
	t64 = t71 * t78 - t81;
	t60 = 0.1e1 / t64 ^ 2;
	t79 = t76 * t73;
	t80 = t74 * t75;
	t63 = t71 * t79 + t80;
	t87 = t60 * t63;
	t86 = t61 * t71;
	t67 = t70 ^ 2;
	t85 = t67 / t71 ^ 2;
	t84 = t70 * t76;
	t66 = 0.1e1 / (t74 ^ 2 * t85 + 0.1e1);
	t83 = t74 * t66;
	t77 = t63 ^ 2 * t60 + 0.1e1;
	t68 = 0.1e1 / t71;
	t59 = 0.1e1 / t64;
	t58 = 0.1e1 / t77;
	t57 = (-0.1e1 - t85) * t83;
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (t67 * t88 + 0.1e1);
	t52 = (-t59 * t73 + t75 * t87) * t58 * t84;
	t51 = (t71 * t54 - (-t74 * t86 - t62 * t70 + (-t62 * t82 - t86) * t57) * t70 * t55) * t76 * t53;
	t1 = [-t68 * t66 * t84, t57, t57, t57, 0; (-t54 * t82 - (t62 * t67 * t68 * t83 + (t66 - 0.1e1) * t70 * t61) * t70 * t88) * t53, t51, t51, t51, 0; ((-t71 * t81 + t78) * t59 - (-t71 * t80 - t79) * t87) * t58, t52, t52, t52, t77 * t58;];
	Ja_rot = t1;
end