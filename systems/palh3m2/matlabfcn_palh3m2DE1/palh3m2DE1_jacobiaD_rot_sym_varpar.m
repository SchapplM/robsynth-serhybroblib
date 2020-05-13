% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m2DE1
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
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh3m2DE1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 02:05
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh3m2DE1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE1_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE1_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2DE1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE1_jacobiaD_rot_sym_varpar: pkin has to be [18x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:34
	% EndTime: 2020-05-07 01:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:34
	% EndTime: 2020-05-07 01:57:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:35
	% EndTime: 2020-05-07 01:57:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:35
	% EndTime: 2020-05-07 01:57:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:36
	% EndTime: 2020-05-07 01:57:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:36
	% EndTime: 2020-05-07 01:57:37
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (2114->46), mult. (3281->122), div. (126->12), fcn. (4544->13), ass. (0->65)
	t168 = sin(pkin(16));
	t169 = cos(pkin(16));
	t172 = sin(pkin(15));
	t175 = cos(pkin(15));
	t161 = t168 * t175 + t169 * t172;
	t162 = -t168 * t172 + t169 * t175;
	t165 = pkin(17) + pkin(18);
	t163 = sin(t165);
	t164 = cos(t165);
	t209 = -t161 * t163 + t164 * t162;
	t156 = 0.1e1 / t209 ^ 2;
	t158 = t161 * t164 + t163 * t162;
	t154 = t158 ^ 2;
	t171 = sin(qJ(1));
	t166 = t171 ^ 2;
	t143 = t154 * t166 * t156 + 0.1e1;
	t174 = cos(qJ(1));
	t167 = t174 ^ 2;
	t200 = 0.1e1 / t143 ^ 2 * t167;
	t208 = t156 * t200;
	t141 = 0.1e1 / t143;
	t207 = (t141 - 0.1e1) * t158;
	t206 = t209 * t174;
	t196 = t158 * t171;
	t140 = atan2(-t196, t209);
	t138 = sin(t140);
	t139 = cos(t140);
	t132 = -t138 * t196 + t139 * t209;
	t129 = 0.1e1 / t132;
	t173 = cos(qJ(4));
	t190 = t174 * t173;
	t170 = sin(qJ(4));
	t193 = t171 * t170;
	t152 = -t190 * t209 + t193;
	t145 = 0.1e1 / t152;
	t155 = 0.1e1 / t209;
	t130 = 0.1e1 / t132 ^ 2;
	t146 = 0.1e1 / t152 ^ 2;
	t184 = t209 * t171;
	t148 = t170 * t184 - t190;
	t133 = t148 * qJD(1) + t152 * qJD(4);
	t192 = t171 * t173;
	t150 = t206 * t170 + t192;
	t144 = t150 ^ 2;
	t137 = t144 * t146 + 0.1e1;
	t199 = t146 * t150;
	t191 = t174 * t170;
	t149 = t173 * t184 + t191;
	t134 = (t191 * t209 + t192) * qJD(4) + t149 * qJD(1);
	t147 = t145 * t146;
	t201 = t134 * t147;
	t204 = (-t133 * t199 - t144 * t201) / t137 ^ 2;
	t203 = t130 * t174;
	t202 = t133 * t146;
	t198 = t149 * t150;
	t197 = t154 * t155;
	t189 = qJD(1) * t171;
	t188 = t155 * t208;
	t124 = (t139 * t141 * t171 * t197 + t138 * t207) * t174;
	t153 = t158 * t154;
	t135 = 0.1e1 / t137;
	t131 = t129 * t130;
	t128 = t154 * t167 * t130 + 0.1e1;
	t123 = qJD(1) * t124;
	t1 = [(t141 * t155 * t158 + 0.2e1 * t153 * t188) * t189, 0, 0, 0; (0.2e1 * (t124 * t203 + t129 * t171) / t128 ^ 2 * (-t123 * t131 * t167 - t189 * t203) * t154 + ((0.2e1 * t124 * t131 * t174 + t130 * t171) * t123 + (-t174 * t129 + ((t124 + (t153 * t208 + t207) * t174 * t138) * t171 - (-0.2e1 * t166 * t154 ^ 2 * t188 + (-t200 + (-t166 + 0.2e1 * t167) * t141) * t197) * t174 * t139) * t130) * qJD(1)) / t128) * t158, 0, 0, 0; 0.2e1 * (-t145 * t148 - t146 * t198) * t204 + (-t149 * t202 + (-t148 * t146 - 0.2e1 * t147 * t198) * t134 + (t149 * t145 + (-t193 * t209 + t190) * t199) * qJD(4) + (t150 * t145 + (t206 * t173 - t193) * t199) * qJD(1)) * t135, 0, 0, -0.2e1 * t204 + 0.2e1 * (-t135 * t202 + (-t135 * t201 - t146 * t204) * t150) * t150;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:37
	% EndTime: 2020-05-07 01:57:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:38
	% EndTime: 2020-05-07 01:57:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:57:39
	% EndTime: 2020-05-07 01:57:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end