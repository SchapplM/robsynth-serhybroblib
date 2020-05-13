% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
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
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh2m1OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh2m1OL_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1OL_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
JaD_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:40
	% EndTime: 2020-05-03 00:53:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 00:53:41
	% EndTime: 2020-05-03 00:53:42
	% DurationCPUTime: 0.84s
	% Computational Cost: add. (7389->96), mult. (5101->211), div. (1026->12), fcn. (5942->9), ass. (0->95)
	t176 = sin(qJ(1));
	t172 = qJ(2) + qJ(3) + qJ(4);
	t169 = sin(t172);
	t165 = t169 ^ 2;
	t170 = cos(t172);
	t167 = 0.1e1 / t170 ^ 2;
	t222 = t165 * t167;
	t195 = 0.1e1 + t222;
	t237 = t176 * t195;
	t173 = t176 ^ 2;
	t161 = t173 * t222 + 0.1e1;
	t159 = 0.1e1 / t161;
	t166 = 0.1e1 / t170;
	t178 = cos(qJ(1));
	t208 = qJD(1) * t178;
	t196 = t169 * t208;
	t171 = qJD(2) + qJD(3) + qJD(4);
	t216 = t171 * t176;
	t199 = t167 * t216;
	t133 = ((-t170 * t216 - t196) * t166 - t165 * t199) * t159;
	t236 = t133 + t216;
	t192 = qJD(1) * t170 + qJD(5);
	t215 = t171 * t178;
	t235 = t169 * t215 + t192 * t176;
	t214 = t176 * t169;
	t158 = atan2(-t214, t170);
	t153 = cos(t158);
	t152 = sin(t158);
	t202 = t152 * t214;
	t143 = t153 * t170 - t202;
	t140 = 0.1e1 / t143;
	t177 = cos(qJ(5));
	t210 = t178 * t177;
	t198 = t170 * t210;
	t175 = sin(qJ(5));
	t213 = t176 * t175;
	t157 = t198 - t213;
	t149 = 0.1e1 / t157;
	t141 = 0.1e1 / t143 ^ 2;
	t150 = 0.1e1 / t157 ^ 2;
	t234 = t159 - 0.1e1;
	t224 = t153 * t169;
	t128 = (-t133 * t176 - t171) * t224 + (-t236 * t170 - t196) * t152;
	t233 = t128 * t140 * t141;
	t138 = -qJD(5) * t198 + t235 * t175 - t177 * t208;
	t211 = t178 * t175;
	t212 = t176 * t177;
	t156 = t170 * t211 + t212;
	t148 = t156 ^ 2;
	t147 = t148 * t150 + 0.1e1;
	t226 = t150 * t156;
	t193 = qJD(5) * t170 + qJD(1);
	t139 = -t235 * t177 - t193 * t211;
	t230 = t139 * t149 * t150;
	t232 = (-t138 * t226 - t148 * t230) / t147 ^ 2;
	t164 = t169 * t165;
	t219 = t166 * t169;
	t188 = t164 * t166 * t167 + t219;
	t220 = t165 * t176;
	t190 = t208 * t220;
	t231 = (t188 * t173 * t171 + t167 * t190) / t161 ^ 2;
	t229 = t141 * t169;
	t228 = t141 * t178;
	t227 = t149 * t175;
	t225 = t152 * t176;
	t223 = t156 * t177;
	t174 = t178 ^ 2;
	t221 = t165 * t174;
	t218 = t169 * t178;
	t217 = t170 * t171;
	t209 = qJD(1) * t176;
	t136 = t141 * t221 + 0.1e1;
	t207 = 0.2e1 * (-t221 * t233 + (t169 * t174 * t217 - t190) * t141) / t136 ^ 2;
	t206 = 0.2e1 * t233;
	t205 = -0.2e1 * t232;
	t204 = t169 * t231;
	t203 = t141 * t218;
	t201 = t156 * t230;
	t200 = t166 * t220;
	t194 = t169 * t207;
	t191 = t153 * t159 * t165 * t166;
	t189 = t195 * t178;
	t187 = t192 * t178;
	t186 = t150 * t223 - t227;
	t155 = -t170 * t212 - t211;
	t154 = -t170 * t213 + t210;
	t145 = 0.1e1 / t147;
	t144 = t159 * t237;
	t134 = 0.1e1 / t136;
	t132 = (t234 * t169 * t152 + t176 * t191) * t178;
	t130 = -t170 * t225 - t224 - (-t152 * t170 - t153 * t214) * t144;
	t129 = 0.2e1 * t231 * t237 + (-qJD(1) * t189 - 0.2e1 * t188 * t216) * t159;
	t126 = t186 * t205 * t218 + (t186 * t170 * t215 + (-t186 * t209 + ((-qJD(5) * t149 - 0.2e1 * t201) * t177 + (-t138 * t177 + (-qJD(5) * t156 + t139) * t175) * t150) * t178) * t169) * t145;
	t125 = (t130 * t229 - t140 * t170) * t178 * t207 + ((-t140 * t209 + (-t130 * t171 - t128) * t228) * t170 + (-t140 * t215 - (-t129 * t153 * t176 + t236 * t152 - (t133 * t225 + t152 * t171 - t153 * t208) * t144) * t203 + (t141 * t209 + t178 * t206) * t130 - ((-t129 - t208) * t152 + ((t144 * t176 - 0.1e1) * t171 + (t144 - t176) * t133) * t153) * t170 * t228) * t169) * t134;
	t1 = [0.2e1 * t178 * t166 * t204 + (-t171 * t189 + t209 * t219) * t159, t129, t129, t129, 0; (t140 * t194 + (-t140 * t217 + (qJD(1) * t132 + t128) * t229) * t134) * t176 + (t141 * t194 * t132 + (-((-0.2e1 * t204 - t217 + (-t133 * t200 + t217) * t159) * t152 + (-0.2e1 * t200 * t231 - t133 * t169 + (t164 * t199 + (t133 + 0.2e1 * t216) * t169) * t159) * t153) * t203 + (-t141 * t217 + t169 * t206) * t132 + (-t140 + ((t173 - t174) * t191 + t234 * t202) * t141) * t169 * qJD(1)) * t134) * t178, t125, t125, t125, 0; 0.2e1 * (-t149 * t154 + t155 * t226) * t232 + (0.2e1 * t155 * t201 - t193 * t149 * t212 + (t171 * t214 - t187) * t227 + (t155 * t138 - t154 * t139 + t187 * t223 - (t169 * t171 * t177 + t193 * t175) * t156 * t176) * t150) * t145, t126, t126, t126, t205 + 0.2e1 * (-t138 * t150 * t145 + (-t145 * t230 - t150 * t232) * t156) * t156;];
	JaD_rot = t1;
end