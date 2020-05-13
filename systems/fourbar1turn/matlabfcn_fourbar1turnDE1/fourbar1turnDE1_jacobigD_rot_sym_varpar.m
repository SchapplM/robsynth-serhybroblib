% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% fourbar1turnDE1
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
% Datum: 2020-04-12 19:28
% Revision: 394980f89398455cc479283a21eae791ed9f69cb (2020-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = fourbar1turnDE1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'fourbar1turnDE1_jacobigD_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'fourbar1turnDE1_jacobigD_rot_sym_varpar: qJD has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1turnDE1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fourbar1turnDE1_jacobigD_rot_sym_varpar: pkin has to be [5x1] (double)');
JgD_rot=NaN(3,2);
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:03
	% EndTime: 2020-04-12 19:28:03
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (2->2), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, qJD(1) * cos(qJ(1)); 0, qJD(1) * sin(qJ(1)); 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:05
	% EndTime: 2020-04-12 19:28:06
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (4924->65), mult. (6850->137), div. (186->9), fcn. (1848->6), ass. (0->71)
	t208 = sin(qJ(2));
	t215 = pkin(1) ^ 2;
	t210 = cos(qJ(2));
	t254 = pkin(1) * t210;
	t240 = -0.2e1 * pkin(2) * t254 + t215;
	t256 = -pkin(3) - pkin(4);
	t198 = (pkin(2) - t256) * (pkin(2) + t256) + t240;
	t255 = pkin(4) - pkin(3);
	t199 = (pkin(2) - t255) * (pkin(2) + t255) + t240;
	t257 = pkin(1) * pkin(2);
	t227 = (-t198 - t199) * t257;
	t191 = t208 * t227;
	t190 = qJD(2) * t191;
	t246 = t198 * t199;
	t216 = sqrt(-t246);
	t194 = 0.1e1 / t216;
	t249 = t194 * t190;
	t248 = t194 * t191;
	t214 = pkin(2) ^ 2;
	t204 = t214 + t240;
	t200 = pkin(3) ^ 2 - pkin(4) ^ 2 + t204;
	t205 = -pkin(2) + t254;
	t244 = t208 * t216;
	t236 = pkin(1) * t244;
	t187 = -t205 * t200 - t236;
	t183 = 0.1e1 / t187;
	t201 = 0.1e1 / t204;
	t184 = 0.1e1 / t187 ^ 2;
	t202 = 0.1e1 / t204 ^ 2;
	t259 = -0.2e1 * t202;
	t258 = -0.2e1 * t205;
	t253 = pkin(3) * t204;
	t243 = t210 * t200;
	t245 = t208 ^ 2 * t215;
	t247 = t194 * t205;
	t179 = -t191 * t247 + 0.2e1 * pkin(2) * t245 + (t243 + t244) * pkin(1);
	t197 = pkin(1) * t208 * t200;
	t188 = -t205 * t216 + t197;
	t213 = 0.1e1 / pkin(3);
	t237 = t208 * t259;
	t226 = t237 * t257;
	t175 = (t179 * t201 + t188 * t226) * t213;
	t252 = t175 * t183;
	t239 = pkin(2) * t258;
	t242 = t216 * t210;
	t176 = (-t208 * t249 + (-t242 + (t200 + t239) * t208) * qJD(2)) * pkin(1);
	t251 = t176 * t183 * t184;
	t250 = t184 * t188;
	t241 = (pkin(1) * t243 + t236) * qJD(2);
	t238 = 0.2e1 * t202;
	t186 = t188 ^ 2;
	t182 = t186 * t184 + 0.1e1;
	t180 = 0.1e1 / t182;
	t235 = t180 * t253;
	t234 = t183 * t253;
	t233 = pkin(2) * qJD(2) * t208;
	t178 = t197 + (-t242 + (t239 - t248) * t208) * pkin(1);
	t174 = (t178 * t201 + t187 * t226) * t213;
	t232 = t174 * t250;
	t231 = qJD(2) * t245;
	t230 = qJD(1) * ((-t232 + t252) * t235 + 0.1e1);
	t228 = pkin(2) * t231;
	t225 = 0.2e1 * t174 * t188 * t253;
	t224 = 0.1e1 / t246 * t190 * t248;
	t222 = 0.8e1 * t201 * t202 * t214 * t231;
	t211 = cos(qJ(1));
	t209 = sin(qJ(1));
	t189 = (t210 * t227 - 0.4e1 * t214 * t245) * qJD(2);
	t177 = -t190 * t247 + 0.2e1 * t228 + t241;
	t171 = (-0.2e1 * t175 * t234 + t184 * t225) * (t177 * t250 - t186 * t251) / t182 ^ 2 + (-t175 * t176 - ((0.4e1 * t228 + t241) * t201 + t187 * t222 + ((-0.2e1 * t210 * t249 + (-t194 * t189 - t224) * t208) * t201 + (t176 * t237 + (t201 * t210 * t258 + (-t178 * t208 - t187 * t210) * t238) * qJD(2)) * pkin(2)) * pkin(1)) * t213 * t188 - t174 * t177) * t184 * t235 + ((0.2e1 * t252 - 0.2e1 * t232) * pkin(1) * pkin(3) * t233 + ((0.6e1 * t215 * t210 * t233 - t189 * t247 - t205 * t224) * t201 + t188 * t222 + ((t177 * pkin(2) * t259 + t201 * t249) * t208 + ((t242 + (-t200 + t248) * t208) * t201 + (-t179 * t208 - t188 * t210) * pkin(2) * t238) * qJD(2)) * pkin(1)) * t213 * t234 + t225 * t251) * t180;
	t1 = [0, t171 * t209 + t211 * t230; 0, -t171 * t211 + t209 * t230; 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-12 19:28:04
	% EndTime: 2020-04-12 19:28:05
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