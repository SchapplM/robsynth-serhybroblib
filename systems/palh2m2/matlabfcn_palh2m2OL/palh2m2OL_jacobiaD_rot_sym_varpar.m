% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
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
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh2m2OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 06:35
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh2m2OL_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'palh2m2OL_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'palh2m2OL_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2OL_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2OL_jacobiaD_rot_sym_varpar: pkin has to be [5x1] (double)');
JaD_rot=NaN(3,6);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:57
	% EndTime: 2020-05-03 06:34:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 06:34:58
	% EndTime: 2020-05-03 06:34:59
	% DurationCPUTime: 0.88s
	% Computational Cost: add. (12423->95), mult. (6392->209), div. (1299->12), fcn. (7429->9), ass. (0->95)
	t169 = sin(qJ(1));
	t166 = t169 ^ 2;
	t165 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
	t162 = sin(t165);
	t157 = t162 ^ 2;
	t163 = cos(t165);
	t159 = 0.1e1 / t163 ^ 2;
	t213 = t157 * t159;
	t154 = t166 * t213 + 0.1e1;
	t156 = t162 * t157;
	t158 = 0.1e1 / t163;
	t164 = qJD(2) + qJD(3) + qJD(4) + qJD(5);
	t210 = t158 * t162;
	t179 = t164 * (t156 * t158 * t159 + t210);
	t171 = cos(qJ(1));
	t199 = qJD(1) * t171;
	t211 = t157 * t169;
	t183 = t199 * t211;
	t222 = (t159 * t183 + t166 * t179) / t154 ^ 2;
	t230 = -0.2e1 * t222;
	t189 = 0.1e1 + t213;
	t229 = t169 * t189;
	t151 = 0.1e1 / t154;
	t190 = t162 * t199;
	t207 = t164 * t169;
	t193 = t159 * t207;
	t126 = ((t207 * t163 + t190) * t158 + t157 * t193) * t151;
	t228 = t126 - t207;
	t185 = qJD(1) * t163 + qJD(6);
	t206 = t164 * t171;
	t227 = t162 * t206 + t169 * t185;
	t205 = t169 * t162;
	t153 = atan2(t205, t163);
	t142 = cos(t153);
	t141 = sin(t153);
	t195 = t141 * t205;
	t136 = t142 * t163 + t195;
	t133 = 0.1e1 / t136;
	t170 = cos(qJ(6));
	t201 = t171 * t170;
	t192 = t163 * t201;
	t168 = sin(qJ(6));
	t204 = t169 * t168;
	t150 = t192 - t204;
	t144 = 0.1e1 / t150;
	t134 = 0.1e1 / t136 ^ 2;
	t145 = 0.1e1 / t150 ^ 2;
	t226 = t151 - 0.1e1;
	t167 = t171 ^ 2;
	t212 = t157 * t167;
	t129 = t134 * t212 + 0.1e1;
	t208 = t163 * t164;
	t217 = t142 * t162;
	t121 = (t126 * t169 - t164) * t217 + (-t228 * t163 + t190) * t141;
	t224 = t121 * t133 * t134;
	t225 = (-t212 * t224 + (t162 * t167 * t208 - t183) * t134) / t129 ^ 2;
	t131 = -qJD(6) * t192 + t227 * t168 - t170 * t199;
	t202 = t171 * t168;
	t203 = t169 * t170;
	t149 = t163 * t202 + t203;
	t143 = t149 ^ 2;
	t140 = t143 * t145 + 0.1e1;
	t215 = t145 * t149;
	t186 = qJD(6) * t163 + qJD(1);
	t132 = -t227 * t170 - t186 * t202;
	t221 = t132 * t144 * t145;
	t223 = (-t131 * t215 - t143 * t221) / t140 ^ 2;
	t220 = t134 * t162;
	t219 = t134 * t171;
	t218 = t141 * t169;
	t216 = t144 * t168;
	t214 = t149 * t170;
	t209 = t162 * t171;
	t200 = qJD(1) * t169;
	t198 = -0.2e1 * t224;
	t197 = -0.2e1 * t223;
	t196 = t134 * t209;
	t194 = t149 * t221;
	t188 = -0.2e1 * t162 * t225;
	t187 = t158 * t230;
	t184 = t142 * t151 * t157 * t158;
	t182 = t189 * t171;
	t181 = t185 * t171;
	t180 = t214 * t145 - t216;
	t148 = -t163 * t203 - t202;
	t147 = -t163 * t204 + t201;
	t138 = 0.1e1 / t140;
	t137 = t151 * t229;
	t127 = 0.1e1 / t129;
	t125 = (-t141 * t162 * t226 + t169 * t184) * t171;
	t123 = t163 * t218 - t217 + (-t141 * t163 + t142 * t205) * t137;
	t122 = t229 * t230 + (qJD(1) * t182 + 0.2e1 * t169 * t179) * t151;
	t119 = t180 * t197 * t209 + (t180 * t163 * t206 + (-t180 * t200 + ((-qJD(6) * t144 - 0.2e1 * t194) * t170 + (-t131 * t170 + (-qJD(6) * t149 + t132) * t168) * t145) * t171) * t162) * t138;
	t118 = 0.2e1 * (-t123 * t220 + t133 * t163) * t171 * t225 + ((t133 * t200 + (t123 * t164 + t121) * t219) * t163 + (t133 * t206 + (t122 * t142 * t169 + t228 * t141 + (-t126 * t218 + t141 * t164 + t142 * t199) * t137) * t196 + (-t134 * t200 + t171 * t198) * t123 + ((-t122 + t199) * t141 + ((t137 * t169 - 0.1e1) * t164 + (-t137 + t169) * t126) * t142) * t163 * t219) * t162) * t127;
	t1 = [t187 * t209 + (t164 * t182 - t200 * t210) * t151, t122, t122, t122, t122, 0; (t133 * t188 + (t133 * t208 + (-qJD(1) * t125 - t121) * t220) * t127) * t169 + (t134 * t188 * t125 + (((0.2e1 * t162 * t222 + t208 + (-t126 * t158 * t211 - t208) * t151) * t141 + (t187 * t211 + t126 * t162 + (t156 * t193 + (-t126 + 0.2e1 * t207) * t162) * t151) * t142) * t196 + (t134 * t208 + t162 * t198) * t125 + (t133 + ((-t166 + t167) * t184 + t226 * t195) * t134) * t162 * qJD(1)) * t127) * t171, t118, t118, t118, t118, 0; 0.2e1 * (-t144 * t147 + t148 * t215) * t223 + (0.2e1 * t148 * t194 - t186 * t144 * t203 + (t164 * t205 - t181) * t216 + (t148 * t131 - t147 * t132 + t181 * t214 - (t162 * t164 * t170 + t168 * t186) * t149 * t169) * t145) * t138, t119, t119, t119, t119, t197 + 0.2e1 * (-t131 * t145 * t138 + (-t138 * t221 - t145 * t223) * t149) * t149;];
	JaD_rot = t1;
end