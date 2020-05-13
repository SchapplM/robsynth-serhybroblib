% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
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
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh3m2TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh3m2TE_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2TE_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_jacobiaD_rot_sym_varpar: pkin has to be [18x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:40
	% EndTime: 2020-05-07 01:41:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:42
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:42
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:43
	% EndTime: 2020-05-07 01:41:43
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (2114->46), mult. (3281->122), div. (126->12), fcn. (4544->13), ass. (0->64)
	t166 = sin(pkin(16));
	t167 = cos(pkin(16));
	t170 = sin(pkin(15));
	t173 = cos(pkin(15));
	t159 = t166 * t173 + t167 * t170;
	t160 = -t166 * t170 + t167 * t173;
	t163 = pkin(17) + pkin(18);
	t161 = sin(t163);
	t162 = cos(t163);
	t206 = -t159 * t161 + t162 * t160;
	t152 = 0.1e1 / t206 ^ 2;
	t157 = t159 * t162 + t160 * t161;
	t155 = t157 ^ 2;
	t169 = sin(qJ(1));
	t164 = t169 ^ 2;
	t141 = t164 * t155 * t152 + 0.1e1;
	t139 = 0.1e1 / t141;
	t205 = (t139 - 0.1e1) * t157;
	t172 = cos(qJ(1));
	t204 = t206 * t172;
	t192 = t169 * t157;
	t138 = atan2(-t192, t206);
	t136 = sin(t138);
	t137 = cos(t138);
	t130 = -t136 * t192 + t137 * t206;
	t127 = 0.1e1 / t130;
	t171 = cos(qJ(4));
	t188 = t172 * t171;
	t168 = sin(qJ(4));
	t191 = t169 * t168;
	t150 = -t188 * t206 + t191;
	t143 = 0.1e1 / t150;
	t151 = 0.1e1 / t206;
	t128 = 0.1e1 / t130 ^ 2;
	t144 = 0.1e1 / t150 ^ 2;
	t182 = t206 * t169;
	t146 = t168 * t182 - t188;
	t131 = t146 * qJD(1) + t150 * qJD(4);
	t190 = t169 * t171;
	t148 = t204 * t168 + t190;
	t142 = t148 ^ 2;
	t135 = t142 * t144 + 0.1e1;
	t197 = t144 * t148;
	t189 = t172 * t168;
	t147 = t171 * t182 + t189;
	t132 = (t189 * t206 + t190) * qJD(4) + t147 * qJD(1);
	t145 = t143 * t144;
	t199 = t132 * t145;
	t202 = (-t131 * t197 - t142 * t199) / t135 ^ 2;
	t201 = t128 * t172;
	t200 = t131 * t144;
	t165 = t172 ^ 2;
	t198 = 0.1e1 / t141 ^ 2 * t165;
	t196 = t147 * t148;
	t195 = t151 * t155;
	t187 = qJD(1) * t169;
	t186 = t157 * t155 * t198;
	t122 = (t137 * t139 * t169 * t195 + t136 * t205) * t172;
	t153 = t151 * t152;
	t133 = 0.1e1 / t135;
	t129 = t127 * t128;
	t126 = t165 * t155 * t128 + 0.1e1;
	t121 = qJD(1) * t122;
	t1 = [(t139 * t151 * t157 + 0.2e1 * t153 * t186) * t187, 0, 0, 0; (0.2e1 * (t122 * t201 + t127 * t169) / t126 ^ 2 * (-t121 * t129 * t165 - t187 * t201) * t155 + ((0.2e1 * t122 * t129 * t172 + t128 * t169) * t121 + (-t172 * t127 + ((t122 + (t152 * t186 + t205) * t172 * t136) * t169 - (-0.2e1 * t153 * t164 * t155 ^ 2 * t198 + (-t198 + (-t164 + 0.2e1 * t165) * t139) * t195) * t172 * t137) * t128) * qJD(1)) / t126) * t157, 0, 0, 0; 0.2e1 * (-t143 * t146 - t144 * t196) * t202 + (-t147 * t200 + (-t146 * t144 - 0.2e1 * t145 * t196) * t132 + (t147 * t143 + (-t191 * t206 + t188) * t197) * qJD(4) + (t148 * t143 + (t204 * t171 - t191) * t197) * qJD(1)) * t133, 0, 0, -0.2e1 * t202 + 0.2e1 * (-t133 * t200 + (-t133 * t199 - t144 * t202) * t148) * t148;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:44
	% EndTime: 2020-05-07 01:41:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end