% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
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
% Datum: 2020-04-12 19:35
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = fourbar1turnDE2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE2_jacobigD_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE2_jacobigD_rot_sym_varpar: qJD has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnDE2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE2_jacobigD_rot_sym_varpar: pkin has to be [5x1] (double)');
JgD_rot=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:42
	% EndTime: 2020-04-12 19:35:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)); 0, qJD(1) * sin(qJ(1)); 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:43
	% EndTime: 2020-04-12 19:35:44
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (4924->65), mult. (6850->137), div. (186->9), fcn. (1848->6), ass. (0->71)
	t159 = sin(qJ(2));
	t166 = pkin(1) ^ 2;
	t161 = cos(qJ(2));
	t205 = pkin(1) * t161;
	t191 = -0.2e1 * pkin(2) * t205 + t166;
	t207 = -pkin(3) - pkin(4);
	t149 = (pkin(2) - t207) * (pkin(2) + t207) + t191;
	t206 = -pkin(3) + pkin(4);
	t150 = (pkin(2) - t206) * (pkin(2) + t206) + t191;
	t208 = pkin(1) * pkin(2);
	t178 = (-t149 - t150) * t208;
	t142 = t159 * t178;
	t141 = qJD(2) * t142;
	t197 = t149 * t150;
	t167 = sqrt(-t197);
	t145 = 0.1e1 / t167;
	t200 = t145 * t141;
	t199 = t145 * t142;
	t165 = pkin(2) ^ 2;
	t155 = t165 + t191;
	t151 = pkin(3) ^ 2 - pkin(4) ^ 2 + t155;
	t156 = -pkin(2) + t205;
	t195 = t159 * t167;
	t187 = pkin(1) * t195;
	t138 = -t156 * t151 - t187;
	t134 = 0.1e1 / t138;
	t152 = 0.1e1 / t155;
	t135 = 0.1e1 / t138 ^ 2;
	t153 = 0.1e1 / t155 ^ 2;
	t210 = -0.2e1 * t153;
	t209 = -0.2e1 * t156;
	t204 = pkin(3) * t155;
	t194 = t161 * t151;
	t196 = t159 ^ 2 * t166;
	t198 = t145 * t156;
	t130 = -t142 * t198 + 0.2e1 * pkin(2) * t196 + (t194 + t195) * pkin(1);
	t148 = pkin(1) * t159 * t151;
	t139 = -t156 * t167 + t148;
	t164 = 0.1e1 / pkin(3);
	t188 = t159 * t210;
	t177 = t188 * t208;
	t126 = (t130 * t152 + t139 * t177) * t164;
	t203 = t126 * t134;
	t190 = pkin(2) * t209;
	t193 = t167 * t161;
	t127 = (-t159 * t200 + (-t193 + (t151 + t190) * t159) * qJD(2)) * pkin(1);
	t202 = t127 * t134 * t135;
	t201 = t135 * t139;
	t192 = (pkin(1) * t194 + t187) * qJD(2);
	t189 = 0.2e1 * t153;
	t137 = t139 ^ 2;
	t133 = t137 * t135 + 0.1e1;
	t131 = 0.1e1 / t133;
	t186 = t131 * t204;
	t185 = t134 * t204;
	t184 = pkin(2) * qJD(2) * t159;
	t129 = t148 + (-t193 + (t190 - t199) * t159) * pkin(1);
	t125 = (t129 * t152 + t138 * t177) * t164;
	t183 = t125 * t201;
	t182 = qJD(2) * t196;
	t181 = qJD(1) * ((-t183 + t203) * t186 + 0.1e1);
	t179 = pkin(2) * t182;
	t176 = 0.2e1 * t125 * t139 * t204;
	t175 = 0.1e1 / t197 * t141 * t199;
	t173 = 0.8e1 * t152 * t153 * t165 * t182;
	t162 = cos(qJ(1));
	t160 = sin(qJ(1));
	t140 = (t161 * t178 - 0.4e1 * t165 * t196) * qJD(2);
	t128 = -t141 * t198 + 0.2e1 * t179 + t192;
	t122 = (-0.2e1 * t126 * t185 + t135 * t176) * (t128 * t201 - t137 * t202) / t133 ^ 2 + (-t126 * t127 - ((0.4e1 * t179 + t192) * t152 + t138 * t173 + ((-0.2e1 * t161 * t200 + (-t145 * t140 - t175) * t159) * t152 + (t127 * t188 + (t152 * t161 * t209 + (-t129 * t159 - t138 * t161) * t189) * qJD(2)) * pkin(2)) * pkin(1)) * t164 * t139 - t125 * t128) * t135 * t186 + ((0.2e1 * t203 - 0.2e1 * t183) * pkin(1) * pkin(3) * t184 + ((0.6e1 * t166 * t161 * t184 - t140 * t198 - t156 * t175) * t152 + t139 * t173 + ((t128 * pkin(2) * t210 + t152 * t200) * t159 + ((t193 + (-t151 + t199) * t159) * t152 + (-t130 * t159 - t139 * t161) * pkin(2) * t189) * qJD(2)) * pkin(1)) * t164 * t185 + t176 * t202) * t131;
	t1 = [0, t122 * t160 + t162 * t181; 0, -t122 * t162 + t160 * t181; 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:35:43
	% EndTime: 2020-04-12 19:35:44
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (4922->63), mult. (6856->142), div. (186->9), fcn. (1846->6), ass. (0->74)
	t179 = pkin(2) ^ 2;
	t180 = pkin(1) ^ 2;
	t175 = cos(qJ(2));
	t220 = pkin(2) * t175;
	t229 = -2 * pkin(1);
	t204 = t220 * t229 + t180;
	t169 = t179 + t204;
	t165 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t169;
	t173 = sin(qJ(2));
	t162 = pkin(2) * t173 * t165;
	t170 = pkin(1) - t220;
	t227 = 2 * pkin(1);
	t200 = t170 * t227;
	t225 = -pkin(3) - pkin(4);
	t163 = (pkin(2) - t225) * (pkin(2) + t225) + t204;
	t224 = -pkin(3) + pkin(4);
	t164 = (pkin(2) - t224) * (pkin(2) + t224) + t204;
	t211 = t163 * t164;
	t181 = sqrt(-t211);
	t206 = t175 * t181;
	t226 = pkin(1) * pkin(2);
	t191 = (-t163 - t164) * t226;
	t156 = t173 * t191;
	t159 = 0.1e1 / t181;
	t213 = t159 * t156;
	t143 = t162 + (-t206 + (t200 - t213) * t173) * pkin(2);
	t207 = t173 * t181;
	t152 = -pkin(2) * t207 + t165 * t170;
	t178 = 0.1e1 / pkin(4);
	t167 = 0.1e1 / t169 ^ 2;
	t209 = t167 * t173;
	t195 = t209 * t226;
	t166 = 0.1e1 / t169;
	t223 = -t166 / 0.2e1;
	t139 = (t143 * t223 + t152 * t195) * t178;
	t149 = 0.1e1 / t152 ^ 2;
	t153 = t170 * t181 + t162;
	t215 = t149 * t153;
	t172 = t173 ^ 2;
	t208 = t172 * t179;
	t210 = t165 * t175;
	t212 = t159 * t170;
	t144 = t156 * t212 + t208 * t227 + (t207 + t210) * pkin(2);
	t222 = t166 / 0.2e1;
	t140 = (t144 * t222 - t153 * t195) * t178;
	t148 = 0.1e1 / t152;
	t217 = t140 * t148;
	t230 = t139 * t215 + t217;
	t155 = qJD(2) * t156;
	t202 = qJD(2) * t179;
	t196 = t172 * t202;
	t192 = pkin(1) * t196;
	t218 = pkin(2) * qJD(2);
	t198 = t173 * t218;
	t205 = t181 * t198 + t210 * t218;
	t142 = t155 * t212 + 0.2e1 * t192 + t205;
	t151 = t153 ^ 2;
	t147 = t149 * t151 + 0.1e1;
	t214 = t159 * t155;
	t141 = (-t173 * t214 + (-t206 + (t165 + t200) * t173) * qJD(2)) * pkin(2);
	t216 = t141 * t148 * t149;
	t228 = (t142 * t215 - t151 * t216) / t147 ^ 2;
	t221 = pkin(1) * t167;
	t219 = pkin(4) * t169;
	t145 = 0.1e1 / t147;
	t199 = t145 * t219;
	t203 = -0.2e1 * qJD(1) * t230 * t199;
	t190 = 0.4e1 / t211 * t155 * t213;
	t188 = t166 * t167 * t180 * t196;
	t176 = cos(qJ(1));
	t174 = sin(qJ(1));
	t154 = (t175 * t191 - 0.4e1 * t180 * t208) * qJD(2);
	t136 = 0.2e1 * t230 * pkin(4) * t145 * t198 * t229 + 0.4e1 * (t217 * t228 + (t145 * t216 + t149 * t228) * t139 * t153) * t219 + 0.2e1 * ((-t139 * t142 + t140 * t141) * t149 + (-((t170 * t190 / 0.4e1 + t154 * t212 + 0.6e1 * t175 * t173 * pkin(1) * t202) * t222 + 0.4e1 * t153 * t188 + ((-t142 * t221 + t214 * t222) * t173 + ((t206 + (-t165 + t213) * t173) * t222 + (-t144 * t173 - t175 * t153) * t221) * qJD(2)) * pkin(2)) * t148 - ((0.4e1 * t192 + t205) * t223 - 0.4e1 * t152 * t188 + ((-0.2e1 * t175 * t214 + (-t190 / 0.4e1 - t159 * t154) * t173) * t223 + (t141 * t209 + (-t166 * t170 * t175 + (t143 * t173 + t152 * t175) * t167) * qJD(2)) * pkin(1)) * pkin(2)) * t215) * t178) * t199;
	t1 = [0, t136 * t174 + t176 * t203; 0, -t136 * t176 + t174 * t203; 0, 0;];
	JgD_rot = t1;
end