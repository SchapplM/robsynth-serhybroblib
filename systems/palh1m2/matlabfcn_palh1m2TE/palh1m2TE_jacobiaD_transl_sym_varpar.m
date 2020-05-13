% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh1m2TE
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh1m2TE_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh1m2TE_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2TE_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_jacobiaD_transl_sym_varpar: pkin has to be [22x1] (double)');
JaD_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:38
	% EndTime: 2020-05-01 20:48:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t18 - r_i_i_C(2) * t16;
	t21 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18 - pkin(15);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t19 * t21) * qJD(1), (-t16 * t26 + t18 * t23) * r_i_i_C(2) + (t16 * t23 + t18 * t26) * r_i_i_C(1), 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t17 * t21) * qJD(1), (t16 * t25 + t18 * t24) * r_i_i_C(2) + (t16 * t24 - t18 * t25) * r_i_i_C(1), 0, 0; 0, -t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (93->25), mult. (246->44), div. (0->0), fcn. (219->6), ass. (0->25)
	t54 = cos(qJ(2));
	t62 = qJD(2) * t54;
	t53 = cos(qJ(3));
	t67 = t53 * t54;
	t50 = sin(qJ(3));
	t51 = sin(qJ(2));
	t68 = t51 * t50;
	t70 = qJD(2) + qJD(3);
	t41 = -qJD(3) * t67 - t53 * t62 + t70 * t68;
	t46 = t54 * t50 + t53 * t51;
	t42 = t70 * t46;
	t66 = -t42 * r_i_i_C(1) + t41 * r_i_i_C(2);
	t58 = pkin(1) * t62 - t66;
	t69 = r_i_i_C(1) * t41;
	t52 = sin(qJ(1));
	t65 = qJD(1) * t52;
	t55 = cos(qJ(1));
	t64 = qJD(1) * t55;
	t63 = qJD(2) * t51;
	t60 = r_i_i_C(1) * qJD(1) * t46;
	t45 = t67 - t68;
	t59 = pkin(1) * t51 - r_i_i_C(1) * t45 + r_i_i_C(2) * t46 - pkin(15);
	t57 = t52 * t60 + t55 * t69 + (t42 * t55 + t45 * t65) * r_i_i_C(2);
	t56 = t52 * t69 - t55 * t60 + (t42 * t52 - t45 * t64) * r_i_i_C(2);
	t1 = [t58 * t52 + (-r_i_i_C(3) * t52 + t59 * t55) * qJD(1), (t54 * t65 + t55 * t63) * pkin(1) + t57, t57, 0; -t58 * t55 + (r_i_i_C(3) * t55 + t59 * t52) * qJD(1), (t52 * t63 - t54 * t64) * pkin(1) + t56, t56, 0; 0, -t58, t66, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (73->29), mult. (178->51), div. (0->0), fcn. (146->12), ass. (0->34)
	t57 = sin(qJ(3));
	t51 = t57 * pkin(5) + pkin(1);
	t62 = cos(qJ(2));
	t75 = t51 * t62;
	t58 = sin(qJ(2));
	t74 = t57 * t58;
	t61 = cos(qJ(3));
	t73 = t58 * t61;
	t72 = t61 * t62;
	t59 = sin(qJ(1));
	t71 = qJD(1) * t59;
	t63 = cos(qJ(1));
	t70 = qJD(1) * t63;
	t69 = qJD(2) * t58;
	t68 = pkin(5) * t72;
	t67 = -t57 * t62 - t73;
	t55 = sin(pkin(20));
	t56 = cos(pkin(20));
	t60 = sin(pkin(18));
	t64 = cos(pkin(18));
	t48 = -t64 * t55 + t60 * t56;
	t49 = t60 * t55 + t64 * t56;
	t54 = pkin(22) + pkin(21);
	t52 = sin(t54);
	t53 = cos(t54);
	t66 = r_i_i_C(1) * (t48 * t52 + t53 * t49) + r_i_i_C(2) * (-t48 * t53 + t52 * t49) + t51 * t58 - pkin(15) - t68;
	t65 = t67 * qJD(3);
	t42 = -qJD(2) * t75 + (-t61 * t69 + t65) * pkin(5);
	t50 = pkin(5) * qJD(3) * t74;
	t47 = t67 * pkin(5);
	t46 = -pkin(5) * t73 - t75;
	t44 = t50 + (-qJD(3) * t72 + (-t72 + t74) * qJD(2)) * pkin(5);
	t43 = t51 * t69 + t50 + (-qJD(2) - qJD(3)) * t68;
	t1 = [-t42 * t59 + (-r_i_i_C(3) * t59 + t66 * t63) * qJD(1), t43 * t63 - t46 * t71, t44 * t63 - t47 * t71, 0; t42 * t63 + (r_i_i_C(3) * t63 + t66 * t59) * qJD(1), t43 * t59 + t46 * t70, t44 * t59 + t47 * t70, 0; 0, t42, (t67 * qJD(2) + t65) * pkin(5), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:40
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (200->47), mult. (442->86), div. (0->0), fcn. (432->14), ass. (0->51)
	t224 = sin(qJ(3));
	t217 = pkin(5) * t224 + pkin(1);
	t230 = cos(qJ(2));
	t251 = t217 * t230;
	t223 = sin(qJ(4));
	t226 = sin(qJ(1));
	t250 = t223 * t226;
	t225 = sin(qJ(2));
	t249 = t224 * t225;
	t229 = cos(qJ(3));
	t248 = t225 * t229;
	t228 = cos(qJ(4));
	t247 = t226 * t228;
	t246 = t229 * t230;
	t245 = qJD(1) * t226;
	t231 = cos(qJ(1));
	t244 = qJD(1) * t231;
	t243 = qJD(2) * t225;
	t242 = pkin(5) * t246;
	t221 = sin(pkin(20));
	t222 = cos(pkin(20));
	t227 = sin(pkin(18));
	t232 = cos(pkin(18));
	t212 = -t221 * t232 + t227 * t222;
	t213 = t227 * t221 + t222 * t232;
	t220 = pkin(22) + pkin(21);
	t218 = sin(t220);
	t219 = cos(t220);
	t241 = t212 * t218 + t213 * t219;
	t240 = -t224 * t230 - t248;
	t239 = t241 * t231;
	t238 = t240 * qJD(3);
	t203 = -t212 * t219 + t213 * t218;
	t214 = -t227 * pkin(9) + pkin(11) * t232;
	t215 = pkin(9) * t232 + t227 * pkin(11);
	t237 = qJD(1) * (-r_i_i_C(3) * t203 - (t214 * t222 + t215 * t221) * t218 + (-t214 * t221 + t215 * t222) * t219 + t217 * t225 - pkin(15) - t242);
	t236 = t228 * t239 - t250;
	t235 = t223 * t231 + t241 * t247;
	t234 = t223 * t239 + t247;
	t233 = t228 * t231 - t241 * t250;
	t204 = -qJD(2) * t251 + (-t229 * t243 + t238) * pkin(5);
	t216 = pkin(5) * qJD(3) * t249;
	t211 = t240 * pkin(5);
	t210 = -pkin(5) * t248 - t251;
	t206 = t216 + (-qJD(3) * t246 + (-t246 + t249) * qJD(2)) * pkin(5);
	t205 = t217 * t243 + t216 + (-qJD(2) - qJD(3)) * t242;
	t202 = t236 * qJD(1) + t233 * qJD(4);
	t201 = t234 * qJD(1) + t235 * qJD(4);
	t200 = t235 * qJD(1) + t234 * qJD(4);
	t199 = t233 * qJD(1) + t236 * qJD(4);
	t1 = [r_i_i_C(1) * t202 - r_i_i_C(2) * t201 - t204 * t226 + t231 * t237, t205 * t231 - t210 * t245, t206 * t231 - t211 * t245, r_i_i_C(1) * t199 - r_i_i_C(2) * t200; r_i_i_C(1) * t200 + r_i_i_C(2) * t199 + t204 * t231 + t226 * t237, t205 * t226 + t210 * t244, t206 * t226 + t211 * t244, r_i_i_C(1) * t201 + r_i_i_C(2) * t202; 0, t204, (t240 * qJD(2) + t238) * pkin(5), (-r_i_i_C(1) * t228 + r_i_i_C(2) * t223) * t203 * qJD(4);];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (71->17), mult. (178->36), div. (0->0), fcn. (200->8), ass. (0->19)
	t59 = sin(pkin(18));
	t62 = cos(pkin(17));
	t69 = sin(pkin(17));
	t70 = cos(pkin(18));
	t54 = t59 * t62 - t70 * t69;
	t55 = t69 * t59 + t62 * t70;
	t57 = sin(qJ(2));
	t60 = cos(qJ(2));
	t65 = t57 * t54 + t55 * t60;
	t46 = t65 * qJD(2);
	t52 = t54 * t60 - t55 * t57;
	t63 = t52 * qJD(2);
	t66 = -r_i_i_C(1) * t46 - r_i_i_C(2) * t63;
	t58 = sin(qJ(1));
	t68 = qJD(1) * t58;
	t61 = cos(qJ(1));
	t67 = qJD(1) * t61;
	t64 = -r_i_i_C(1) * t52 + r_i_i_C(2) * t65 + pkin(14);
	t1 = [-t66 * t58 + (-r_i_i_C(3) * t58 + t64 * t61) * qJD(1), (t46 * t61 + t52 * t68) * r_i_i_C(2) + (-t61 * t63 + t65 * t68) * r_i_i_C(1), 0, 0; t66 * t61 + (r_i_i_C(3) * t61 + t64 * t58) * qJD(1), (t46 * t58 - t52 * t67) * r_i_i_C(2) + (-t58 * t63 - t65 * t67) * r_i_i_C(1), 0, 0; 0, t66, 0, 0;];
	JaD_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->12), mult. (50->24), div. (0->0), fcn. (39->8), ass. (0->13)
	t31 = cos(qJ(2));
	t37 = qJD(1) * t31;
	t28 = sin(qJ(2));
	t36 = qJD(2) * t28;
	t35 = pkin(1) * qJD(2) * t31;
	t26 = sin(pkin(22));
	t27 = cos(pkin(22));
	t30 = sin(pkin(18));
	t33 = cos(pkin(18));
	t34 = r_i_i_C(1) * (t30 * t26 + t27 * t33) + r_i_i_C(2) * (t33 * t26 - t30 * t27) + pkin(1) * t28 - pkin(15);
	t32 = cos(qJ(1));
	t29 = sin(qJ(1));
	t1 = [t29 * t35 + (-r_i_i_C(3) * t29 + t34 * t32) * qJD(1), (t29 * t37 + t32 * t36) * pkin(1), 0, 0; -t32 * t35 + (r_i_i_C(3) * t32 + t34 * t29) * qJD(1), (t29 * t36 - t32 * t37) * pkin(1), 0, 0; 0, -t35, 0, 0;];
	JaD_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (172->31), mult. (432->60), div. (0->0), fcn. (468->8), ass. (0->30)
	t104 = sin(qJ(3));
	t107 = cos(qJ(3));
	t102 = sin(pkin(19));
	t103 = cos(pkin(19));
	t105 = sin(qJ(2));
	t108 = cos(qJ(2));
	t136 = -t102 * t105 + t103 * t108;
	t98 = t102 * t108 + t103 * t105;
	t82 = t104 * t98 - t107 * t136;
	t91 = t98 * qJD(2);
	t92 = t136 * qJD(2);
	t134 = qJD(3) * t82 + t91 * t104 - t92 * t107;
	t117 = t104 * t136 + t98 * t107;
	t113 = -qJD(3) * t117 - t104 * t92 - t91 * t107;
	t139 = r_i_i_C(1) * t113;
	t135 = t134 * r_i_i_C(2) + t139;
	t137 = -t104 * t102 + t107 * t103;
	t95 = t102 * t107 + t103 * t104;
	t138 = t105 * t95 - t108 * t137;
	t106 = sin(qJ(1));
	t115 = -t105 * t137 - t108 * t95;
	t109 = cos(qJ(1));
	t119 = qJD(1) * t109;
	t93 = t95 * qJD(3);
	t94 = t137 * qJD(3);
	t133 = (t138 * t119 + t106 * (-qJD(2) * t115 + t105 * t94 + t108 * t93)) * r_i_i_C(2) + (t106 * t134 - t117 * t119) * r_i_i_C(1);
	t120 = qJD(1) * t106;
	t132 = (-t109 * t113 - t82 * t120) * r_i_i_C(2) + (t109 * t134 + t117 * t120) * r_i_i_C(1);
	t118 = t82 * r_i_i_C(1) - pkin(15);
	t1 = [-t135 * t106 + (-r_i_i_C(3) * t106 + (r_i_i_C(2) * t117 + t118) * t109) * qJD(1), t132, t132, 0; (t139 + (qJD(2) * t138 + t105 * t93 - t94 * t108) * r_i_i_C(2)) * t109 + (t109 * r_i_i_C(3) + (-r_i_i_C(2) * t115 + t118) * t106) * qJD(1), t133, t133, 0; 0, t135, t135, 0;];
	JaD_transl = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (102->27), mult. (350->49), div. (0->0), fcn. (290->8), ass. (0->32)
	t44 = sin(pkin(19));
	t45 = cos(pkin(19));
	t46 = sin(qJ(3));
	t49 = cos(qJ(3));
	t40 = t49 * t44 + t46 * t45;
	t38 = t40 * qJD(3);
	t47 = sin(qJ(2));
	t50 = cos(qJ(2));
	t60 = t46 * t44 - t49 * t45;
	t55 = t60 * t50;
	t72 = t40 * t47 + t55;
	t77 = t72 * qJD(2) + t38 * t47;
	t69 = t40 * t50;
	t75 = t60 * t47;
	t74 = -t69 + t75;
	t36 = t74 * pkin(2);
	t76 = qJD(2) * t74;
	t33 = (qJD(3) * t75 - t38 * t50 + t76) * pkin(2);
	t62 = r_i_i_C(1) * t50 - r_i_i_C(2) * t47;
	t57 = t62 * qJD(2);
	t73 = t33 - t57;
	t48 = sin(qJ(1));
	t65 = qJD(1) * t48;
	t51 = cos(qJ(1));
	t64 = qJD(1) * t51;
	t61 = r_i_i_C(1) * t47 + r_i_i_C(2) * t50;
	t59 = t72 * pkin(2) - pkin(15) + t61;
	t58 = t62 - t36;
	t54 = t61 * qJD(2) + (qJD(3) * t55 + t77) * pkin(2);
	t39 = t60 * qJD(3);
	t34 = (t39 * t50 + t77) * pkin(2);
	t1 = [-t73 * t48 + (-r_i_i_C(3) * t48 + t59 * t51) * qJD(1), t54 * t51 + t58 * t65, t34 * t51 - t36 * t65, 0; t73 * t51 + (r_i_i_C(3) * t51 + t59 * t48) * qJD(1), t54 * t48 - t58 * t64, t34 * t48 + t36 * t64, 0; 0, -t57 + (-qJD(3) * t69 + t39 * t47 + t76) * pkin(2), t33, 0;];
	JaD_transl = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (85->29), mult. (130->44), div. (0->0), fcn. (93->12), ass. (0->30)
	t42 = qJD(2) + qJD(3);
	t43 = qJ(3) + qJ(2);
	t41 = cos(t43);
	t66 = r_i_i_C(2) * t41;
	t40 = sin(t43);
	t68 = r_i_i_C(1) * t40;
	t57 = t66 + t68;
	t55 = t57 * t42;
	t50 = cos(qJ(2));
	t69 = pkin(1) * t50;
	t70 = qJD(2) * t69 + t55;
	t67 = r_i_i_C(2) * t40;
	t65 = t41 * t42;
	t48 = sin(qJ(1));
	t64 = qJD(1) * t48;
	t51 = cos(qJ(1));
	t63 = qJD(1) * t51;
	t47 = sin(qJ(2));
	t62 = qJD(2) * t47;
	t61 = r_i_i_C(1) * t65;
	t60 = t42 * t67;
	t58 = qJD(1) * t66;
	t44 = sin(pkin(22));
	t46 = cos(pkin(22));
	t49 = sin(pkin(18));
	t52 = cos(pkin(18));
	t56 = -r_i_i_C(1) * t41 + pkin(1) * t47 - pkin(15) - (-(t49 * t44 + t46 * t52) * cos(pkin(21)) + (t52 * t44 - t49 * t46) * sin(pkin(21))) * pkin(4) + t67;
	t54 = t48 * t58 + t64 * t68 + (t60 - t61) * t51;
	t36 = t48 * t60;
	t1 = [t70 * t48 + (-r_i_i_C(3) * t48 + t56 * t51) * qJD(1), (t50 * t64 + t51 * t62) * pkin(1) + t54, t54, 0; -t70 * t51 + (r_i_i_C(3) * t51 + t56 * t48) * qJD(1), t36 + (pkin(1) * t62 - t61) * t48 + (-t57 - t69) * t63, -t51 * t58 + t36 + (-t40 * t63 - t48 * t65) * r_i_i_C(1), 0; 0, -t70, -t55, 0;];
	JaD_transl = t1;
end