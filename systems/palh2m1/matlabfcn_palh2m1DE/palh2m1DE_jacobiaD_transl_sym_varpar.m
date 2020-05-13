% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh2m1DE
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh2m1DE_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh2m1DE_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1DE_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
JaD_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->13), mult. (60->29), div. (0->0), fcn. (38->4), ass. (0->12)
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t19 = cos(qJ(1));
	t25 = qJD(1) * t19;
	t24 = qJD(2) * t17;
	t23 = qJD(2) * t19;
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(1);
	t20 = t22 * qJD(2);
	t1 = [t17 * t20 + (r_i_i_C(3) * t17 + t21 * t19) * qJD(1), (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0; -t22 * t23 + (-r_i_i_C(3) * t19 + t21 * t17) * qJD(1), (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0; 0, t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (77->25), mult. (110->40), div. (0->0), fcn. (71->6), ass. (0->27)
	t41 = qJD(2) + qJD(3);
	t42 = qJ(2) + qJ(3);
	t40 = cos(t42);
	t60 = r_i_i_C(2) * t40;
	t39 = sin(t42);
	t61 = r_i_i_C(1) * t39;
	t49 = t60 + t61;
	t43 = sin(qJ(2));
	t62 = pkin(2) * t43;
	t50 = qJD(2) * t62;
	t63 = t49 * t41 + t50;
	t59 = t39 * t41;
	t58 = t40 * t41;
	t57 = r_i_i_C(1) * t59 + r_i_i_C(2) * t58;
	t44 = sin(qJ(1));
	t56 = qJD(1) * t44;
	t46 = cos(qJ(1));
	t55 = qJD(1) * t46;
	t45 = cos(qJ(2));
	t54 = qJD(2) * t45;
	t53 = r_i_i_C(1) * t58;
	t52 = r_i_i_C(2) * t59;
	t51 = qJD(1) * t60;
	t48 = -t45 * pkin(2) - r_i_i_C(1) * t40 + r_i_i_C(2) * t39 - pkin(1);
	t47 = t44 * t51 + t56 * t61 + (t52 - t53) * t46;
	t32 = t44 * t52;
	t1 = [t63 * t44 + (r_i_i_C(3) * t44 + t48 * t46) * qJD(1), (t43 * t56 - t46 * t54) * pkin(2) + t47, t47, 0; -t63 * t46 + (-r_i_i_C(3) * t46 + t48 * t44) * qJD(1), t32 + (-pkin(2) * t54 - t53) * t44 + (-t49 - t62) * t55, -t46 * t51 + t32 + (-t39 * t55 - t44 * t58) * r_i_i_C(1), 0; 0, t50 + t57, t57, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (43->16), mult. (62->30), div. (0->0), fcn. (39->6), ass. (0->19)
	t22 = qJ(2) + qJ(3);
	t19 = sin(t22);
	t33 = pkin(3) * t19;
	t20 = cos(t22);
	t25 = cos(qJ(2));
	t32 = -t25 * pkin(2) - pkin(3) * t20 - pkin(1) - r_i_i_C(1);
	t21 = qJD(2) + qJD(3);
	t31 = t20 * t21;
	t30 = pkin(2) * qJD(2);
	t24 = sin(qJ(1));
	t29 = qJD(1) * t24;
	t26 = cos(qJ(1));
	t28 = qJD(1) * t26;
	t27 = t21 * t33;
	t23 = sin(qJ(2));
	t18 = -t23 * pkin(2) - t33;
	t16 = -pkin(3) * t31 - t25 * t30;
	t15 = t23 * t30 + t27;
	t1 = [t24 * t15 + (r_i_i_C(3) * t24 + t32 * t26) * qJD(1), t26 * t16 - t18 * t29, (t19 * t29 - t26 * t31) * pkin(3), 0; -t26 * t15 + (-r_i_i_C(3) * t26 + t32 * t24) * qJD(1), t24 * t16 + t18 * t28, (-t19 * t28 - t24 * t31) * pkin(3), 0; 0, t15, t27, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (71->24), mult. (138->41), div. (0->0), fcn. (107->8), ass. (0->27)
	t100 = qJD(1) + qJD(4);
	t80 = qJ(2) + qJ(3);
	t77 = sin(t80);
	t99 = pkin(3) * t77;
	t79 = qJD(2) + qJD(3);
	t83 = sin(qJ(1));
	t98 = t79 * t83;
	t86 = cos(qJ(1));
	t97 = t79 * t86;
	t81 = sin(qJ(4));
	t84 = cos(qJ(4));
	t72 = t100 * (t81 * t83 - t84 * t86);
	t73 = t100 * (-t81 * t86 - t83 * t84);
	t96 = t72 * r_i_i_C(1) - t73 * r_i_i_C(2);
	t95 = t73 * r_i_i_C(1) + t72 * r_i_i_C(2);
	t94 = qJD(1) * t83;
	t93 = qJD(1) * t86;
	t85 = cos(qJ(2));
	t92 = qJD(2) * t85;
	t82 = sin(qJ(2));
	t91 = qJD(2) * t82 * pkin(2);
	t78 = cos(t80);
	t88 = -pkin(3) * t78 * t97 + t94 * t99;
	t87 = (-t77 * t93 - t78 * t98) * pkin(3);
	t76 = t85 * pkin(2) + pkin(1) + pkin(4);
	t75 = t79 * t99;
	t1 = [t83 * t91 - t76 * t93 + (t77 * t98 - t78 * t93) * pkin(3) + t96, (t82 * t94 - t86 * t92) * pkin(2) + t88, t88, t96; -t86 * t91 - t76 * t94 + (-t77 * t97 - t78 * t94) * pkin(3) + t95, t87 + (-t82 * t93 - t83 * t92) * pkin(2), t87, t95; 0, t75 + t91, t75, 0;];
	JaD_transl = t1;
end