% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m2DE2
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
%   Wie in palh3m2DE2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% Ja_rot [3x4]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = palh3m2DE2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_jacobia_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_jacobia_rot_sym_varpar: pkin has to be [18x1] (double)');
Ja_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:22
	% EndTime: 2020-05-07 02:13:22
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:23
	% EndTime: 2020-05-07 02:13:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:24
	% EndTime: 2020-05-07 02:13:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:30
	% EndTime: 2020-05-07 02:13:31
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (10602->23), mult. (18186->50), div. (211->9), fcn. (27825->20), ass. (0->43)
	t167 = pkin(17) + pkin(18);
	t165 = sin(t167);
	t166 = cos(t167);
	t172 = sin(qJ(2));
	t177 = cos(qJ(2));
	t168 = sin(pkin(16));
	t169 = cos(pkin(16));
	t174 = sin(pkin(15));
	t179 = cos(pkin(15));
	t163 = t168 * t179 + t169 * t174;
	t164 = -t168 * t174 + t169 * t179;
	t171 = sin(qJ(3));
	t176 = cos(qJ(3));
	t184 = t171 * t163 - t164 * t176;
	t186 = -t163 * t176 - t171 * t164;
	t185 = t172 * t184 + t177 * t186;
	t208 = -t172 * t186 + t184 * t177;
	t148 = qJ(2) + qJ(3) + atan2(t165 * t208 + t166 * t185, t185 * t165 - t208 * t166);
	t146 = sin(t148);
	t143 = t146 ^ 2;
	t147 = cos(t148);
	t173 = sin(qJ(1));
	t136 = 0.1e1 / (t173 ^ 2 * t143 / t147 ^ 2 + 0.1e1);
	t216 = t136 / t147;
	t192 = t173 * t146;
	t137 = atan2(t192, t147);
	t134 = sin(t137);
	t135 = cos(t137);
	t130 = t134 * t192 + t135 * t147;
	t178 = cos(qJ(1));
	t204 = 0.1e1 / t130 ^ 2 * t178 ^ 2;
	t170 = sin(qJ(4));
	t191 = t173 * t170;
	t175 = cos(qJ(4));
	t190 = t173 * t175;
	t189 = t178 * t170;
	t188 = t178 * t175;
	t142 = -t147 * t188 + t191;
	t139 = 0.1e1 / t142 ^ 2;
	t140 = t147 * t189 + t190;
	t187 = t140 ^ 2 * t139 + 0.1e1;
	t133 = 0.1e1 / t187;
	t1 = [t146 * t178 * t216, 0, 0, 0; (0.1e1 / t130 * t192 + (t135 * t143 * t173 * t216 + (-t136 + 0.1e1) * t146 * t134) * t146 * t204) / (t143 * t204 + 0.1e1), 0, 0, 0; ((t147 * t191 - t188) / t142 + (t147 * t190 + t189) * t139 * t140) * t133, 0, 0, t187 * t133;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobia_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:52
	% EndTime: 2020-05-07 02:13:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobia_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 02:13:54
	% EndTime: 2020-05-07 02:13:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
end