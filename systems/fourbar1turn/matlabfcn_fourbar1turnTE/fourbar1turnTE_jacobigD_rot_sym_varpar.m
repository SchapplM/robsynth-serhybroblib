% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
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
% qJD [2x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4,l5]';
% 
% Output:
% JgD_rot [3x2]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-12 19:20
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = fourbar1turnTE_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnTE_jacobigD_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnTE_jacobigD_rot_sym_varpar: qJD has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnTE_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnTE_jacobigD_rot_sym_varpar: pkin has to be [5x1] (double)');
JgD_rot=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:25
	% EndTime: 2020-04-12 19:20:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)); 0, qJD(1) * sin(qJ(1)); 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:26
	% EndTime: 2020-04-12 19:20:27
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (4924->65), mult. (6850->137), div. (186->9), fcn. (1848->6), ass. (0->71)
	t189 = sin(qJ(2));
	t196 = pkin(1) ^ 2;
	t191 = cos(qJ(2));
	t235 = pkin(1) * t191;
	t221 = -0.2e1 * pkin(2) * t235 + t196;
	t237 = -pkin(3) - pkin(4);
	t179 = (pkin(2) - t237) * (pkin(2) + t237) + t221;
	t236 = -pkin(3) + pkin(4);
	t180 = (pkin(2) - t236) * (pkin(2) + t236) + t221;
	t238 = pkin(1) * pkin(2);
	t208 = (-t179 - t180) * t238;
	t172 = t189 * t208;
	t171 = qJD(2) * t172;
	t227 = t179 * t180;
	t197 = sqrt(-t227);
	t175 = 0.1e1 / t197;
	t230 = t175 * t171;
	t229 = t175 * t172;
	t195 = pkin(2) ^ 2;
	t185 = t195 + t221;
	t181 = pkin(3) ^ 2 - pkin(4) ^ 2 + t185;
	t186 = -pkin(2) + t235;
	t224 = t189 * t197;
	t217 = pkin(1) * t224;
	t168 = -t181 * t186 - t217;
	t164 = 0.1e1 / t168;
	t182 = 0.1e1 / t185;
	t165 = 0.1e1 / t168 ^ 2;
	t183 = 0.1e1 / t185 ^ 2;
	t240 = -0.2e1 * t183;
	t239 = -0.2e1 * t186;
	t234 = pkin(3) * t185;
	t225 = t189 ^ 2 * t196;
	t226 = t181 * t191;
	t228 = t175 * t186;
	t160 = -t172 * t228 + 0.2e1 * pkin(2) * t225 + (t224 + t226) * pkin(1);
	t178 = pkin(1) * t189 * t181;
	t169 = -t186 * t197 + t178;
	t194 = 0.1e1 / pkin(3);
	t218 = t189 * t240;
	t207 = t218 * t238;
	t156 = (t160 * t182 + t169 * t207) * t194;
	t233 = t156 * t164;
	t220 = pkin(2) * t239;
	t223 = t197 * t191;
	t157 = (-t189 * t230 + (-t223 + (t181 + t220) * t189) * qJD(2)) * pkin(1);
	t232 = t157 * t164 * t165;
	t231 = t165 * t169;
	t222 = (pkin(1) * t226 + t217) * qJD(2);
	t219 = 0.2e1 * t183;
	t167 = t169 ^ 2;
	t163 = t165 * t167 + 0.1e1;
	t161 = 0.1e1 / t163;
	t216 = t161 * t234;
	t215 = t164 * t234;
	t214 = pkin(2) * qJD(2) * t189;
	t159 = t178 + (-t223 + (t220 - t229) * t189) * pkin(1);
	t155 = (t159 * t182 + t168 * t207) * t194;
	t213 = t155 * t231;
	t212 = qJD(2) * t225;
	t211 = qJD(1) * ((-t213 + t233) * t216 + 0.1e1);
	t209 = pkin(2) * t212;
	t206 = 0.2e1 * t155 * t169 * t234;
	t205 = 0.1e1 / t227 * t171 * t229;
	t203 = 0.8e1 * t182 * t183 * t195 * t212;
	t192 = cos(qJ(1));
	t190 = sin(qJ(1));
	t170 = (t191 * t208 - 0.4e1 * t195 * t225) * qJD(2);
	t158 = -t171 * t228 + 0.2e1 * t209 + t222;
	t152 = (-0.2e1 * t156 * t215 + t165 * t206) * (t158 * t231 - t167 * t232) / t163 ^ 2 + (-t156 * t157 - ((0.4e1 * t209 + t222) * t182 + t168 * t203 + ((-0.2e1 * t191 * t230 + (-t175 * t170 - t205) * t189) * t182 + (t157 * t218 + (t182 * t191 * t239 + (-t159 * t189 - t168 * t191) * t219) * qJD(2)) * pkin(2)) * pkin(1)) * t194 * t169 - t155 * t158) * t165 * t216 + ((0.2e1 * t233 - 0.2e1 * t213) * pkin(1) * pkin(3) * t214 + ((0.6e1 * t196 * t191 * t214 - t170 * t228 - t186 * t205) * t182 + t169 * t203 + ((t158 * pkin(2) * t240 + t182 * t230) * t189 + ((t223 + (-t181 + t229) * t189) * t182 + (-t160 * t189 - t169 * t191) * pkin(2) * t219) * qJD(2)) * pkin(1)) * t194 * t215 + t206 * t232) * t161;
	t1 = [0, t152 * t190 + t192 * t211; 0, -t152 * t192 + t190 * t211; 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:20:26
	% EndTime: 2020-04-12 19:20:27
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (4922->63), mult. (6856->142), div. (186->9), fcn. (1846->6), ass. (0->74)
	t157 = pkin(2) ^ 2;
	t158 = pkin(1) ^ 2;
	t153 = cos(qJ(2));
	t198 = pkin(2) * t153;
	t207 = -2 * pkin(1);
	t182 = t198 * t207 + t158;
	t147 = t157 + t182;
	t143 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t147;
	t151 = sin(qJ(2));
	t140 = pkin(2) * t151 * t143;
	t148 = pkin(1) - t198;
	t205 = 2 * pkin(1);
	t178 = t148 * t205;
	t203 = -pkin(3) - pkin(4);
	t141 = (pkin(2) - t203) * (pkin(2) + t203) + t182;
	t202 = -pkin(3) + pkin(4);
	t142 = (pkin(2) - t202) * (pkin(2) + t202) + t182;
	t189 = t141 * t142;
	t159 = sqrt(-t189);
	t184 = t153 * t159;
	t204 = pkin(1) * pkin(2);
	t169 = (-t141 - t142) * t204;
	t134 = t151 * t169;
	t137 = 0.1e1 / t159;
	t191 = t137 * t134;
	t121 = t140 + (-t184 + (t178 - t191) * t151) * pkin(2);
	t185 = t151 * t159;
	t130 = -pkin(2) * t185 + t143 * t148;
	t156 = 0.1e1 / pkin(4);
	t145 = 0.1e1 / t147 ^ 2;
	t187 = t145 * t151;
	t173 = t187 * t204;
	t144 = 0.1e1 / t147;
	t201 = -t144 / 0.2e1;
	t117 = (t121 * t201 + t130 * t173) * t156;
	t127 = 0.1e1 / t130 ^ 2;
	t131 = t148 * t159 + t140;
	t193 = t127 * t131;
	t150 = t151 ^ 2;
	t186 = t150 * t157;
	t188 = t143 * t153;
	t190 = t137 * t148;
	t122 = t134 * t190 + t186 * t205 + (t185 + t188) * pkin(2);
	t200 = t144 / 0.2e1;
	t118 = (t122 * t200 - t131 * t173) * t156;
	t126 = 0.1e1 / t130;
	t195 = t118 * t126;
	t208 = t117 * t193 + t195;
	t133 = qJD(2) * t134;
	t180 = qJD(2) * t157;
	t174 = t150 * t180;
	t170 = pkin(1) * t174;
	t196 = pkin(2) * qJD(2);
	t176 = t151 * t196;
	t183 = t159 * t176 + t188 * t196;
	t120 = t133 * t190 + 0.2e1 * t170 + t183;
	t129 = t131 ^ 2;
	t125 = t127 * t129 + 0.1e1;
	t192 = t137 * t133;
	t119 = (-t151 * t192 + (-t184 + (t143 + t178) * t151) * qJD(2)) * pkin(2);
	t194 = t119 * t126 * t127;
	t206 = (t120 * t193 - t129 * t194) / t125 ^ 2;
	t199 = pkin(1) * t145;
	t197 = pkin(4) * t147;
	t123 = 0.1e1 / t125;
	t177 = t123 * t197;
	t181 = -0.2e1 * qJD(1) * t208 * t177;
	t168 = 0.4e1 / t189 * t133 * t191;
	t166 = t144 * t145 * t158 * t174;
	t154 = cos(qJ(1));
	t152 = sin(qJ(1));
	t132 = (t153 * t169 - 0.4e1 * t158 * t186) * qJD(2);
	t114 = 0.2e1 * t208 * pkin(4) * t123 * t176 * t207 + 0.4e1 * (t195 * t206 + (t123 * t194 + t127 * t206) * t117 * t131) * t197 + 0.2e1 * ((-t117 * t120 + t118 * t119) * t127 + (-((t148 * t168 / 0.4e1 + t132 * t190 + 0.6e1 * t153 * t151 * pkin(1) * t180) * t200 + 0.4e1 * t131 * t166 + ((-t120 * t199 + t192 * t200) * t151 + ((t184 + (-t143 + t191) * t151) * t200 + (-t122 * t151 - t153 * t131) * t199) * qJD(2)) * pkin(2)) * t126 - ((0.4e1 * t170 + t183) * t201 - 0.4e1 * t130 * t166 + ((-0.2e1 * t153 * t192 + (-t168 / 0.4e1 - t137 * t132) * t151) * t201 + (t119 * t187 + (-t144 * t148 * t153 + (t121 * t151 + t130 * t153) * t145) * qJD(2)) * pkin(1)) * pkin(2)) * t193) * t156) * t177;
	t1 = [0, t114 * t152 + t154 * t181; 0, -t114 * t154 + t152 * t181; 0, 0;];
	JgD_rot = t1;
end