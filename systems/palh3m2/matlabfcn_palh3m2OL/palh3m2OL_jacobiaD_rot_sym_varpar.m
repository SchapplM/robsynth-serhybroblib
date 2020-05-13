% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% palh3m2OL
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% qJD [10x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in palh3m2OL_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% JaD_rot [3x10]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:44
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = palh3m2OL_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(10,1),uint8(0),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m2OL_jacobiaD_rot_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [10 1]), ...
  'palh3m2OL_jacobiaD_rot_sym_varpar: qJD has to be [10x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2OL_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m2OL_jacobiaD_rot_sym_varpar: pkin has to be [16x1] (double)');
JaD_rot=NaN(3,10);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:52
	% EndTime: 2020-05-07 04:44:53
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (7008->94), mult. (5101->206), div. (1026->12), fcn. (5942->9), ass. (0->94)
	t171 = sin(qJ(1));
	t168 = t171 ^ 2;
	t167 = qJ(2) + qJ(3) + qJ(4);
	t164 = sin(t167);
	t160 = t164 ^ 2;
	t165 = cos(t167);
	t162 = 0.1e1 / t165 ^ 2;
	t214 = t160 * t162;
	t158 = t168 * t214 + 0.1e1;
	t159 = t164 * t160;
	t161 = 0.1e1 / t165;
	t166 = qJD(2) + qJD(3) + qJD(4);
	t211 = t161 * t164;
	t182 = t166 * (t159 * t161 * t162 + t211);
	t173 = cos(qJ(1));
	t200 = qJD(1) * t173;
	t212 = t160 * t171;
	t186 = t200 * t212;
	t223 = (t162 * t186 + t168 * t182) / t158 ^ 2;
	t230 = -0.2e1 * t223;
	t193 = 0.1e1 + t214;
	t229 = t171 * t193;
	t155 = 0.1e1 / t158;
	t194 = t164 * t200;
	t208 = t166 * t171;
	t195 = t162 * t208;
	t129 = ((t165 * t208 + t194) * t161 + t160 * t195) * t155;
	t228 = t129 - t208;
	t206 = t171 * t164;
	t157 = atan2(t206, t165);
	t154 = cos(t157);
	t153 = sin(t157);
	t196 = t153 * t206;
	t139 = t154 * t165 + t196;
	t136 = 0.1e1 / t139;
	t172 = cos(qJ(5));
	t202 = t173 * t172;
	t170 = sin(qJ(5));
	t205 = t171 * t170;
	t152 = -t165 * t202 + t205;
	t145 = 0.1e1 / t152;
	t137 = 0.1e1 / t139 ^ 2;
	t146 = 0.1e1 / t152 ^ 2;
	t227 = t155 - 0.1e1;
	t169 = t173 ^ 2;
	t213 = t160 * t169;
	t132 = t137 * t213 + 0.1e1;
	t209 = t165 * t166;
	t215 = t154 * t164;
	t124 = (t129 * t171 - t166) * t215 + (-t228 * t165 + t194) * t153;
	t225 = t124 * t136 * t137;
	t226 = (-t213 * t225 + (t164 * t169 * t209 - t186) * t137) / t132 ^ 2;
	t188 = qJD(1) * t165 + qJD(5);
	t207 = t166 * t173;
	t181 = t164 * t207 + t188 * t171;
	t189 = qJD(5) * t165 + qJD(1);
	t134 = t181 * t170 - t189 * t202;
	t203 = t173 * t170;
	t204 = t171 * t172;
	t150 = t165 * t203 + t204;
	t144 = t150 ^ 2;
	t143 = t144 * t146 + 0.1e1;
	t218 = t146 * t150;
	t135 = t181 * t172 + t189 * t203;
	t222 = t135 * t145 * t146;
	t224 = (-t134 * t218 - t144 * t222) / t143 ^ 2;
	t221 = t137 * t164;
	t220 = t137 * t173;
	t219 = t145 * t170;
	t217 = t150 * t172;
	t216 = t153 * t171;
	t210 = t164 * t173;
	t201 = qJD(1) * t171;
	t199 = -0.2e1 * t225;
	t198 = -0.2e1 * t224;
	t197 = t137 * t210;
	t192 = -0.2e1 * t164 * t226;
	t191 = t161 * t230;
	t190 = -0.2e1 * t150 * t222;
	t187 = t154 * t155 * t160 * t161;
	t185 = t193 * t173;
	t184 = t188 * t173;
	t183 = t146 * t217 + t219;
	t149 = t165 * t204 + t203;
	t148 = t165 * t205 - t202;
	t141 = 0.1e1 / t143;
	t140 = t155 * t229;
	t130 = 0.1e1 / t132;
	t128 = (-t227 * t164 * t153 + t171 * t187) * t173;
	t126 = t165 * t216 - t215 + (-t153 * t165 + t154 * t206) * t140;
	t125 = t229 * t230 + (qJD(1) * t185 + 0.2e1 * t171 * t182) * t155;
	t122 = t183 * t198 * t210 + (t183 * t165 * t207 + (-t183 * t201 + ((qJD(5) * t145 + t190) * t172 + (-t134 * t172 + (-qJD(5) * t150 - t135) * t170) * t146) * t173) * t164) * t141;
	t121 = 0.2e1 * (-t126 * t221 + t136 * t165) * t173 * t226 + ((t136 * t201 + (t126 * t166 + t124) * t220) * t165 + (t136 * t207 + (t125 * t154 * t171 + t228 * t153 + (-t129 * t216 + t153 * t166 + t154 * t200) * t140) * t197 + (-t137 * t201 + t173 * t199) * t126 + ((-t125 + t200) * t153 + ((t140 * t171 - 0.1e1) * t166 + (-t140 + t171) * t129) * t154) * t165 * t220) * t164) * t130;
	t1 = [t191 * t210 + (t166 * t185 - t201 * t211) * t155, t125, t125, t125, 0, 0, 0, 0, 0, 0; (t136 * t192 + (t136 * t209 + (-qJD(1) * t128 - t124) * t221) * t130) * t171 + (t137 * t192 * t128 + (((0.2e1 * t164 * t223 + t209 + (-t129 * t161 * t212 - t209) * t155) * t153 + (t191 * t212 + t129 * t164 + (t159 * t195 + (-t129 + 0.2e1 * t208) * t164) * t155) * t154) * t197 + (t137 * t209 + t164 * t199) * t128 + (t136 + ((-t168 + t169) * t187 + t227 * t196) * t137) * t164 * qJD(1)) * t130) * t173, t121, t121, t121, 0, 0, 0, 0, 0, 0; 0.2e1 * (-t145 * t148 - t149 * t218) * t224 + (t149 * t190 + t189 * t145 * t204 + (-t166 * t206 + t184) * t219 + (-t149 * t134 - t148 * t135 + t184 * t217 + (-t164 * t166 * t172 - t189 * t170) * t150 * t171) * t146) * t141, t122, t122, t122, t198 + 0.2e1 * (-t134 * t146 * t141 + (-t141 * t222 - t146 * t224) * t150) * t150, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 04:44:51
	% EndTime: 2020-05-07 04:44:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
end