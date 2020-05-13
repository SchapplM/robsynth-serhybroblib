% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh3m2TE
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
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh3m2TE_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh3m2TE_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2TE_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_jacobiaD_transl_sym_varpar: pkin has to be [18x1] (double)');
JaD_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:40
	% EndTime: 2020-05-07 01:41:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
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
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->14), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(12);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0; 0, -t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:42
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (93->22), mult. (246->38), div. (0->0), fcn. (219->6), ass. (0->23)
	t58 = sin(qJ(3));
	t59 = sin(qJ(2));
	t61 = cos(qJ(3));
	t62 = cos(qJ(2));
	t55 = t58 * t62 + t61 * t59;
	t79 = qJD(2) + qJD(3);
	t50 = t79 * t55;
	t54 = t58 * t59 - t61 * t62;
	t51 = t79 * t54;
	t72 = t50 * r_i_i_C(1) - t51 * r_i_i_C(2);
	t78 = pkin(1) * t59;
	t64 = qJD(2) * t78 - t72;
	t66 = -r_i_i_C(1) * t55 + r_i_i_C(2) * t54;
	t80 = -r_i_i_C(1) * t51 - r_i_i_C(2) * t50;
	t63 = cos(qJ(1));
	t73 = t80 * t63;
	t60 = sin(qJ(1));
	t71 = qJD(1) * t60;
	t70 = qJD(1) * t63;
	t69 = qJD(2) * t62;
	t67 = t80 * t60 - t66 * t70;
	t65 = -pkin(1) * t62 - r_i_i_C(1) * t54 - r_i_i_C(2) * t55 - pkin(12);
	t1 = [t64 * t60 + (-r_i_i_C(3) * t60 + t65 * t63) * qJD(1), -pkin(1) * t63 * t69 + (t66 + t78) * t71 + t73, t66 * t71 + t73, 0; -t64 * t63 + (r_i_i_C(3) * t63 + t65 * t60) * qJD(1), (-t59 * t70 - t60 * t69) * pkin(1) + t67, t67, 0; 0, -t64, t72, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:42
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->30), mult. (178->50), div. (0->0), fcn. (146->12), ass. (0->35)
	t67 = cos(qJ(3));
	t81 = pkin(4) * t67;
	t68 = cos(qJ(2));
	t80 = pkin(4) * t68;
	t63 = sin(qJ(3));
	t64 = sin(qJ(2));
	t79 = t63 * t64;
	t65 = sin(qJ(1));
	t78 = qJD(1) * t65;
	t69 = cos(qJ(1));
	t77 = qJD(1) * t69;
	t57 = pkin(1) - t81;
	t76 = qJD(2) * t57;
	t56 = t63 * t80;
	t75 = pkin(4) * t79;
	t74 = t64 * t81;
	t72 = qJD(3) * t80;
	t73 = qJD(2) * t56 + qJD(3) * t74 + t63 * t72;
	t61 = sin(pkin(16));
	t62 = cos(pkin(16));
	t66 = sin(pkin(15));
	t70 = cos(pkin(15));
	t50 = t61 * t70 + t62 * t66;
	t51 = -t61 * t66 + t62 * t70;
	t60 = pkin(17) + pkin(18);
	t58 = sin(t60);
	t59 = cos(t60);
	t71 = -r_i_i_C(1) * (t58 * t50 - t51 * t59) + r_i_i_C(2) * (t50 * t59 + t51 * t58) - t57 * t68 - pkin(12) - t75;
	t55 = t67 * t72;
	t49 = t56 + t74;
	t48 = -t57 * t64 + t56;
	t46 = t55 + (-qJD(3) * t79 + (t67 * t68 - t79) * qJD(2)) * pkin(4);
	t45 = -t64 * t76 + t73;
	t44 = -t68 * t76 + t55 + (-qJD(2) - qJD(3)) * t75;
	t1 = [-t45 * t65 + (-r_i_i_C(3) * t65 + t71 * t69) * qJD(1), t44 * t69 - t48 * t78, t46 * t69 - t49 * t78, 0; t45 * t69 + (r_i_i_C(3) * t69 + t71 * t65) * qJD(1), t44 * t65 + t48 * t77, t46 * t65 + t49 * t77, 0; 0, t45, qJD(2) * t74 + t73, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:43
	% EndTime: 2020-05-07 01:41:44
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (200->48), mult. (442->85), div. (0->0), fcn. (432->14), ass. (0->52)
	t241 = cos(qJ(3));
	t263 = pkin(4) * t241;
	t242 = cos(qJ(2));
	t262 = pkin(4) * t242;
	t235 = sin(qJ(4));
	t238 = sin(qJ(1));
	t261 = t235 * t238;
	t236 = sin(qJ(3));
	t237 = sin(qJ(2));
	t260 = t236 * t237;
	t240 = cos(qJ(4));
	t259 = t238 * t240;
	t258 = qJD(1) * t238;
	t243 = cos(qJ(1));
	t257 = qJD(1) * t243;
	t229 = pkin(1) - t263;
	t256 = qJD(2) * t229;
	t255 = pkin(4) * t260;
	t228 = t236 * t262;
	t254 = t237 * t263;
	t253 = qJD(3) * t262;
	t252 = qJD(2) * t228 + qJD(3) * t254 + t236 * t253;
	t233 = sin(pkin(16));
	t234 = cos(pkin(16));
	t239 = sin(pkin(15));
	t244 = cos(pkin(15));
	t220 = t233 * t244 + t234 * t239;
	t221 = -t233 * t239 + t234 * t244;
	t232 = pkin(17) + pkin(18);
	t230 = sin(t232);
	t231 = cos(t232);
	t251 = t220 * t230 - t221 * t231;
	t250 = t251 * t243;
	t211 = t220 * t231 + t221 * t230;
	t222 = t239 * pkin(8) + pkin(10) * t244;
	t223 = pkin(8) * t244 - t239 * pkin(10);
	t249 = qJD(1) * (-r_i_i_C(3) * t211 - (t222 * t234 + t223 * t233) * t230 + (-t222 * t233 + t223 * t234) * t231 - t229 * t242 - pkin(12) - t255);
	t248 = -t240 * t250 - t261;
	t247 = t235 * t243 - t251 * t259;
	t246 = -t235 * t250 + t259;
	t245 = t240 * t243 + t251 * t261;
	t227 = t241 * t253;
	t219 = t228 + t254;
	t218 = -t229 * t237 + t228;
	t214 = t227 + (-qJD(3) * t260 + (t241 * t242 - t260) * qJD(2)) * pkin(4);
	t213 = -t237 * t256 + t252;
	t212 = -t242 * t256 + t227 + (-qJD(2) - qJD(3)) * t255;
	t210 = t248 * qJD(1) + t245 * qJD(4);
	t209 = t246 * qJD(1) + t247 * qJD(4);
	t208 = t247 * qJD(1) + t246 * qJD(4);
	t207 = t245 * qJD(1) + t248 * qJD(4);
	t1 = [r_i_i_C(1) * t210 - r_i_i_C(2) * t209 - t213 * t238 + t243 * t249, t212 * t243 - t218 * t258, t214 * t243 - t219 * t258, r_i_i_C(1) * t207 - r_i_i_C(2) * t208; t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t213 * t243 + t238 * t249, t212 * t238 + t218 * t257, t214 * t238 + t219 * t257, r_i_i_C(1) * t209 + r_i_i_C(2) * t210; 0, t213, qJD(2) * t254 + t252, (-r_i_i_C(1) * t240 + r_i_i_C(2) * t235) * t211 * qJD(4);];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:44
	% EndTime: 2020-05-07 01:41:44
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (71->17), mult. (178->36), div. (0->0), fcn. (200->8), ass. (0->19)
	t59 = sin(pkin(15));
	t62 = cos(pkin(14));
	t69 = sin(pkin(14));
	t70 = cos(pkin(15));
	t54 = t62 * t59 - t70 * t69;
	t55 = t69 * t59 + t62 * t70;
	t57 = sin(qJ(2));
	t60 = cos(qJ(2));
	t65 = t60 * t54 + t57 * t55;
	t48 = t65 * qJD(2);
	t53 = t57 * t54 - t55 * t60;
	t63 = t53 * qJD(2);
	t66 = -r_i_i_C(1) * t48 + r_i_i_C(2) * t63;
	t58 = sin(qJ(1));
	t68 = qJD(1) * t58;
	t61 = cos(qJ(1));
	t67 = qJD(1) * t61;
	t64 = r_i_i_C(1) * t53 + r_i_i_C(2) * t65 + pkin(6);
	t1 = [-t66 * t58 + (-r_i_i_C(3) * t58 + t64 * t61) * qJD(1), (t48 * t61 - t53 * t68) * r_i_i_C(2) + (t61 * t63 + t65 * t68) * r_i_i_C(1), 0, 0; t66 * t61 + (r_i_i_C(3) * t61 + t64 * t58) * qJD(1), (t48 * t58 + t53 * t67) * r_i_i_C(2) + (t58 * t63 - t65 * t67) * r_i_i_C(1), 0, 0; 0, t66, 0, 0;];
	JaD_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (17->12), mult. (50->24), div. (0->0), fcn. (39->8), ass. (0->13)
	t27 = sin(qJ(2));
	t36 = qJD(1) * t27;
	t30 = cos(qJ(2));
	t35 = qJD(2) * t30;
	t34 = pkin(1) * qJD(2) * t27;
	t25 = sin(pkin(18));
	t26 = cos(pkin(18));
	t29 = sin(pkin(15));
	t32 = cos(pkin(15));
	t33 = r_i_i_C(1) * (-t25 * t29 + t26 * t32) + r_i_i_C(2) * (t25 * t32 + t26 * t29) - pkin(1) * t30 - pkin(12);
	t31 = cos(qJ(1));
	t28 = sin(qJ(1));
	t1 = [t28 * t34 + (-r_i_i_C(3) * t28 + t33 * t31) * qJD(1), (t28 * t36 - t31 * t35) * pkin(1), 0, 0; -t31 * t34 + (r_i_i_C(3) * t31 + t33 * t28) * qJD(1), (-t28 * t35 - t31 * t36) * pkin(1), 0, 0; 0, -t34, 0, 0;];
	JaD_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (85->29), mult. (130->47), div. (0->0), fcn. (93->12), ass. (0->31)
	t51 = qJD(2) + qJD(3);
	t52 = qJ(3) + qJ(2);
	t50 = cos(t52);
	t75 = r_i_i_C(2) * t50;
	t49 = sin(t52);
	t76 = r_i_i_C(1) * t49;
	t64 = t75 + t76;
	t56 = sin(qJ(2));
	t77 = pkin(1) * t56;
	t66 = qJD(2) * t77;
	t78 = t64 * t51 - t66;
	t74 = t49 * t51;
	t73 = t50 * t51;
	t72 = r_i_i_C(1) * t74 + r_i_i_C(2) * t73;
	t57 = sin(qJ(1));
	t71 = qJD(1) * t57;
	t60 = cos(qJ(1));
	t70 = qJD(1) * t60;
	t59 = cos(qJ(2));
	t69 = qJD(2) * t59;
	t68 = r_i_i_C(1) * t73;
	t67 = r_i_i_C(2) * t74;
	t65 = qJD(1) * t76;
	t53 = sin(pkin(18));
	t55 = cos(pkin(18));
	t58 = sin(pkin(15));
	t61 = cos(pkin(15));
	t63 = r_i_i_C(1) * t50 - r_i_i_C(2) * t49 - pkin(1) * t59 - pkin(12) - (-(-t53 * t58 + t55 * t61) * cos(pkin(17)) + (t61 * t53 + t58 * t55) * sin(pkin(17))) * pkin(3);
	t62 = t60 * t65 + t70 * t75 + (-t67 + t68) * t57;
	t44 = t60 * t68;
	t1 = [-t78 * t57 + (-r_i_i_C(3) * t57 + t63 * t60) * qJD(1), t44 + (-pkin(1) * t69 - t67) * t60 + (-t64 + t77) * t71, -t57 * t65 + t44 + (-t50 * t71 - t60 * t74) * r_i_i_C(2), 0; t78 * t60 + (r_i_i_C(3) * t60 + t63 * t57) * qJD(1), (-t56 * t70 - t57 * t69) * pkin(1) + t62, t62, 0; 0, -t66 + t72, t72, 0;];
	JaD_transl = t1;
end