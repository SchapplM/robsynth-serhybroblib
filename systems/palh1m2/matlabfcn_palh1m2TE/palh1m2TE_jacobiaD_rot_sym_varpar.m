% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh1m2TE
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
%   Wie in palh1m2TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh1m2TE_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2TE_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_jacobiaD_rot_sym_varpar: pkin has to be [22x1] (double)');
JaD_rot=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:38
	% EndTime: 2020-05-01 20:48:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:40
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (2114->45), mult. (3281->121), div. (126->12), fcn. (4544->13), ass. (0->64)
	t159 = sin(qJ(4));
	t160 = sin(qJ(1));
	t162 = cos(qJ(4));
	t181 = t160 * t162;
	t157 = sin(pkin(20));
	t158 = cos(pkin(20));
	t161 = sin(pkin(18));
	t164 = cos(pkin(18));
	t150 = -t164 * t157 + t161 * t158;
	t151 = t161 * t157 + t164 * t158;
	t154 = pkin(22) + pkin(21);
	t152 = sin(t154);
	t153 = cos(t154);
	t148 = t150 * t152 + t153 * t151;
	t163 = cos(qJ(1));
	t197 = t148 * t163;
	t140 = t159 * t197 + t181;
	t182 = t160 * t159;
	t198 = t162 * t197 - t182;
	t144 = 0.1e1 / t148 ^ 2;
	t149 = -t150 * t153 + t152 * t151;
	t147 = t149 ^ 2;
	t155 = t160 ^ 2;
	t133 = t155 * t147 * t144 + 0.1e1;
	t131 = 0.1e1 / t133;
	t196 = (t131 - 0.1e1) * t149;
	t183 = t160 * t149;
	t130 = atan2(-t183, t148);
	t128 = sin(t130);
	t129 = cos(t130);
	t122 = -t128 * t183 + t129 * t148;
	t119 = 0.1e1 / t122;
	t135 = 0.1e1 / t198;
	t143 = 0.1e1 / t148;
	t120 = 0.1e1 / t122 ^ 2;
	t136 = 0.1e1 / t198 ^ 2;
	t173 = t148 * t160;
	t179 = t163 * t162;
	t138 = t159 * t173 - t179;
	t123 = t138 * qJD(1) - qJD(4) * t198;
	t134 = t140 ^ 2;
	t127 = t134 * t136 + 0.1e1;
	t188 = t136 * t140;
	t180 = t163 * t159;
	t139 = t162 * t173 + t180;
	t124 = (t148 * t180 + t181) * qJD(4) + t139 * qJD(1);
	t137 = t135 * t136;
	t190 = t124 * t137;
	t193 = (-t123 * t188 + t134 * t190) / t127 ^ 2;
	t192 = t120 * t163;
	t191 = t123 * t136;
	t156 = t163 ^ 2;
	t189 = 0.1e1 / t133 ^ 2 * t156;
	t187 = t139 * t140;
	t186 = t143 * t147;
	t178 = qJD(1) * t160;
	t177 = t149 * t147 * t189;
	t114 = (t129 * t131 * t160 * t186 + t128 * t196) * t163;
	t145 = t143 * t144;
	t125 = 0.1e1 / t127;
	t121 = t119 * t120;
	t118 = t156 * t147 * t120 + 0.1e1;
	t113 = qJD(1) * t114;
	t1 = [(t131 * t143 * t149 + 0.2e1 * t145 * t177) * t178, 0, 0, 0; (0.2e1 * (t114 * t192 + t119 * t160) / t118 ^ 2 * (-t113 * t121 * t156 - t178 * t192) * t147 + ((0.2e1 * t114 * t121 * t163 + t120 * t160) * t113 + (-t163 * t119 + ((t114 + (t144 * t177 + t196) * t163 * t128) * t160 - (-0.2e1 * t145 * t155 * t147 ^ 2 * t189 + (-t189 + (-t155 + 0.2e1 * t156) * t131) * t186) * t163 * t129) * t120) * qJD(1)) / t118) * t149, 0, 0, 0; 0.2e1 * (t135 * t138 - t136 * t187) * t193 + (-t139 * t191 + (-t138 * t136 + 0.2e1 * t137 * t187) * t124 + (-t139 * t135 + (-t148 * t182 + t179) * t188) * qJD(4) + (-t140 * t135 + t198 * t188) * qJD(1)) * t125, 0, 0, -0.2e1 * t193 + 0.2e1 * (-t125 * t191 + (t125 * t190 - t136 * t193) * t140) * t140;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
end