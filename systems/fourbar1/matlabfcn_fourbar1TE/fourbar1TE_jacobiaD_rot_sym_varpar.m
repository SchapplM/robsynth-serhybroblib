% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in fourbar1TE_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% JaD_rot [3x1]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = fourbar1TE_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_jacobiaD_rot_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_jacobiaD_rot_sym_varpar: qJD has to be [1x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'fourbar1TE_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_jacobiaD_rot_sym_varpar: pkin has to be [4x1] (double)');
JaD_rot=NaN(3,1);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0; 0; 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:42
	% EndTime: 2020-04-24 19:49:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0; 0; 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:43
	% EndTime: 2020-04-24 19:49:44
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (2178->62), mult. (3042->125), div. (79->9), fcn. (815->4), ass. (0->64)
	t173 = sin(qJ(1));
	t178 = pkin(1) ^ 2;
	t174 = cos(qJ(1));
	t211 = pkin(2) * t174;
	t218 = -2 * pkin(1);
	t198 = t211 * t218 + t178;
	t216 = -pkin(3) - pkin(4);
	t163 = (pkin(2) - t216) * (pkin(2) + t216) + t198;
	t215 = -pkin(3) + pkin(4);
	t164 = (pkin(2) - t215) * (pkin(2) + t215) + t198;
	t217 = pkin(1) * pkin(2);
	t189 = (-t163 - t164) * t217;
	t157 = t173 * t189;
	t156 = qJD(1) * t157;
	t177 = pkin(2) ^ 2;
	t169 = t177 + t198;
	t165 = pkin(3) ^ 2 - pkin(4) ^ 2 + t169;
	t203 = t165 * t174;
	t161 = pkin(2) * qJD(1) * t203;
	t201 = t173 ^ 2 * t177;
	t192 = 0.2e1 * pkin(1) * t201;
	t204 = t163 * t164;
	t179 = sqrt(-t204);
	t200 = t173 * t179;
	t196 = pkin(2) * t200;
	t159 = 0.1e1 / t179;
	t170 = -pkin(1) + t211;
	t205 = t159 * t170;
	t143 = t156 * t205 + t161 + (t192 - t196) * qJD(1);
	t153 = -t165 * t170 + t196;
	t150 = 0.1e1 / t153 ^ 2;
	t162 = pkin(2) * t173 * t165;
	t154 = t170 * t179 + t162;
	t152 = t154 ^ 2;
	t148 = t150 * t152 + 0.1e1;
	t208 = t150 * t154;
	t197 = t170 * t218;
	t199 = t174 * t179;
	t207 = t159 * t156;
	t142 = (t173 * t207 + (t199 + (t165 + t197) * t173) * qJD(1)) * pkin(2);
	t149 = 0.1e1 / t153;
	t209 = t142 * t149 * t150;
	t219 = (t143 * t208 - t152 * t209) / t148 ^ 2;
	t206 = t159 * t157;
	t166 = 0.1e1 / t169;
	t167 = 0.1e1 / t169 ^ 2;
	t214 = -t166 / 0.2e1;
	t213 = t166 / 0.2e1;
	t212 = pkin(1) * t167;
	t210 = pkin(3) * t169;
	t202 = t167 * t173;
	t195 = t149 * t210;
	t194 = pkin(1) * qJD(1) * t173;
	t193 = qJD(1) * t201;
	t191 = t202 * t217;
	t146 = 0.1e1 / t148;
	t190 = t146 * t150 * t210;
	t187 = t166 * t167 * t178 * t193;
	t186 = 0.1e1 / t204 * t156 * t206;
	t185 = 0.2e1 * pkin(2) * pkin(3) * t146 * t194;
	t155 = (t174 * t189 - 0.4e1 * t178 * t201) * qJD(1);
	t145 = t157 * t205 + t192 + (-t200 + t203) * pkin(2);
	t144 = t162 + (t199 + (t197 + t206) * t173) * pkin(2);
	t1 = [0; 0; 0.2e1 * (((0.4e1 * pkin(1) * t193 + t161) * t213 + 0.4e1 * t153 * t187 + ((0.2e1 * t174 * t207 + (-qJD(1) * t179 + t159 * t155 + t186) * t173) * t213 + (-t142 * t202 + (-t166 * t170 * t174 + (-t144 * t173 - t153 * t174) * t167) * qJD(1)) * pkin(1)) * pkin(2)) * t154 * t190 + ((0.6e1 * t177 * t174 * t194 + t155 * t205 + t170 * t186) * t214 - 0.4e1 * t154 * t187 + ((t143 * t212 + t207 * t213) * t173 + ((-t199 + (-t165 - t206) * t173) * t214 + (t145 * t173 + t174 * t154) * t212) * qJD(1)) * pkin(2)) * t146 * t195 + (-t142 * t190 + t149 * t185 - 0.2e1 * t195 * t219) * (t145 * t214 + t154 * t191) + (t143 * t190 + t185 * t208 - 0.2e1 * (t146 * t209 + t150 * t219) * t154 * t210) * (t144 * t213 - t153 * t191)) / pkin(3);];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-04-24 19:49:43
	% EndTime: 2020-04-24 19:49:43
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (2178->59), mult. (3042->126), div. (79->9), fcn. (815->4), ass. (0->63)
	t171 = sin(qJ(1));
	t176 = pkin(1) ^ 2;
	t172 = cos(qJ(1));
	t207 = pkin(2) * t172;
	t214 = -2 * pkin(1);
	t196 = t207 * t214 + t176;
	t212 = -pkin(3) - pkin(4);
	t161 = (pkin(2) - t212) * (pkin(2) + t212) + t196;
	t211 = -pkin(3) + pkin(4);
	t162 = (pkin(2) - t211) * (pkin(2) + t211) + t196;
	t213 = pkin(1) * pkin(2);
	t188 = (-t161 - t162) * t213;
	t157 = t171 * t188;
	t156 = qJD(1) * t157;
	t175 = pkin(2) ^ 2;
	t167 = t175 + t196;
	t163 = -pkin(3) ^ 2 + pkin(4) ^ 2 + t167;
	t201 = t161 * t162;
	t177 = sqrt(-t201);
	t198 = t171 * t177;
	t185 = -t172 * t163 - t198;
	t199 = t171 ^ 2 * t175;
	t183 = t185 * pkin(2) + t199 * t214;
	t159 = 0.1e1 / t177;
	t168 = -pkin(1) + t207;
	t202 = t159 * t168;
	t143 = t183 * qJD(1) + t156 * t202;
	t153 = pkin(2) * t198 + t163 * t168;
	t150 = 0.1e1 / t153 ^ 2;
	t154 = -pkin(2) * t171 * t163 + t168 * t177;
	t152 = t154 ^ 2;
	t148 = t150 * t152 + 0.1e1;
	t204 = t150 * t154;
	t195 = 0.2e1 * t168 * pkin(1);
	t197 = t172 * t177;
	t203 = t159 * t156;
	t142 = (t171 * t203 + (t197 + (-t163 + t195) * t171) * qJD(1)) * pkin(2);
	t149 = 0.1e1 / t153;
	t205 = t142 * t149 * t150;
	t217 = (t143 * t204 - t152 * t205) / t148 ^ 2;
	t216 = t159 * t157;
	t215 = t172 * t159;
	t164 = 0.1e1 / t167;
	t165 = 0.1e1 / t167 ^ 2;
	t210 = -t164 / 0.2e1;
	t209 = t164 / 0.2e1;
	t208 = pkin(1) * t165;
	t206 = pkin(4) * t167;
	t200 = t165 * t171;
	t194 = t149 * t206;
	t193 = pkin(1) * qJD(1) * t171;
	t192 = t200 * t213;
	t191 = -0.4e1 * t176 * t199;
	t190 = t163 - t216;
	t146 = 0.1e1 / t148;
	t189 = t146 * t150 * t206;
	t186 = 0.1e1 / t201 * t156 * t216;
	t184 = 0.2e1 * pkin(2) * pkin(4) * t146 * t193;
	t166 = t164 * t165;
	t155 = (t172 * t188 + t191) * qJD(1);
	t145 = t157 * t202 + t183;
	t144 = (t197 + (-t190 + t195) * t171) * pkin(2);
	t1 = [0; 0; 0.2e1 * (((0.4e1 * t153 * t166 * t176 + t164 * t214) * qJD(1) * t199 + ((t156 * t215 + (t157 * t215 + t185) * qJD(1) + (t159 * t155 + t186) * t171) * t209 + (-t142 * t200 + (t164 * t168 * t172 + (-t144 * t171 - t153 * t172) * t165) * qJD(1)) * pkin(1)) * pkin(2)) * t154 * t189 + ((-0.6e1 * t175 * t172 * t193 + t155 * t202 + t168 * t186) * t210 + t154 * t166 * qJD(1) * t191 + ((t143 * t208 + t203 * t209) * t171 + ((t190 * t171 - t197) * t210 + (t145 * t171 + t172 * t154) * t208) * qJD(1)) * pkin(2)) * t146 * t194 + (-t142 * t189 + t149 * t184 - 0.2e1 * t194 * t217) * (t145 * t210 + t154 * t192) + (t143 * t189 + t184 * t204 - 0.2e1 * (t205 * t146 + t150 * t217) * t154 * t206) * (t144 * t209 - t153 * t192)) / pkin(4);];
	JaD_rot = t1;
end