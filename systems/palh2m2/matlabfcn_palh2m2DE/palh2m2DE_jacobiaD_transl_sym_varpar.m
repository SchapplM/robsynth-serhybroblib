% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% palh2m2DE
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[A2Off,A3Off,A4Off,L1,L2]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 01:06
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = palh2m2DE_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m2DE_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m2DE_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'palh2m2DE_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m2DE_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'palh2m2DE_jacobiaD_transl_sym_varpar: pkin has to be [5x1] (double)');
JaD_transl=NaN(3,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
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
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
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
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18;
	t21 = -r_i_i_C(1) * t18 + r_i_i_C(2) * t16 - pkin(1);
	t20 = t22 * qJD(2);
	t1 = [t22 * t24 + (-r_i_i_C(3) * t17 + t21 * t19) * qJD(1), (t16 * t23 + t18 * t26) * r_i_i_C(2) + (t16 * t26 - t18 * t23) * r_i_i_C(1), 0, 0; -t19 * t20 + (r_i_i_C(3) * t19 + t21 * t17) * qJD(1), (t16 * t24 - t18 * t25) * r_i_i_C(2) + (-t16 * t25 - t18 * t24) * r_i_i_C(1), 0, 0; 0, -t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->9), mult. (34->18), div. (0->0), fcn. (21->4), ass. (0->9)
	t13 = cos(qJ(2));
	t18 = -t13 * pkin(4) - pkin(1) - r_i_i_C(1);
	t11 = sin(qJ(2));
	t17 = qJD(1) * t11;
	t16 = qJD(2) * t13;
	t15 = qJD(2) * t11 * pkin(4);
	t14 = cos(qJ(1));
	t12 = sin(qJ(1));
	t1 = [t12 * t15 + (-r_i_i_C(3) * t12 + t18 * t14) * qJD(1), (t12 * t17 - t14 * t16) * pkin(4), 0, 0; -t14 * t15 + (r_i_i_C(3) * t14 + t18 * t12) * qJD(1), (-t12 * t16 - t14 * t17) * pkin(4), 0, 0; 0, -t15, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (26->20), mult. (82->38), div. (0->0), fcn. (53->6), ass. (0->16)
	t21 = sin(qJ(3));
	t24 = cos(qJ(3));
	t27 = (r_i_i_C(1) * t21 + r_i_i_C(2) * t24) * qJD(3);
	t22 = sin(qJ(2));
	t30 = qJD(2) * t22 * pkin(4);
	t36 = -t27 - t30;
	t23 = sin(qJ(1));
	t35 = qJD(1) * t23;
	t26 = cos(qJ(1));
	t34 = qJD(1) * t26;
	t25 = cos(qJ(2));
	t33 = qJD(2) * t25;
	t32 = qJD(3) * t23;
	t31 = qJD(3) * t26;
	t28 = -t25 * pkin(4) - r_i_i_C(1) * t24 + r_i_i_C(2) * t21 - pkin(1) - pkin(2);
	t1 = [-t36 * t23 + (-r_i_i_C(3) * t23 + t26 * t28) * qJD(1), (t22 * t35 - t26 * t33) * pkin(4), (t21 * t31 + t24 * t35) * r_i_i_C(2) + (t21 * t35 - t24 * t31) * r_i_i_C(1), 0; t36 * t26 + (r_i_i_C(3) * t26 + t23 * t28) * qJD(1), (-t22 * t34 - t23 * t33) * pkin(4), (t21 * t32 - t24 * t34) * r_i_i_C(2) + (-t21 * t34 - t24 * t32) * r_i_i_C(1), 0; 0, -t30, -t27, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (20->15), mult. (54->28), div. (0->0), fcn. (34->6), ass. (0->15)
	t18 = cos(qJ(3));
	t19 = cos(qJ(2));
	t27 = -t19 * pkin(4) - t18 * pkin(5) - pkin(1) - pkin(2) - r_i_i_C(1);
	t17 = sin(qJ(1));
	t26 = qJD(1) * t17;
	t20 = cos(qJ(1));
	t25 = qJD(1) * t20;
	t24 = qJD(2) * t19;
	t23 = qJD(3) * t18;
	t16 = sin(qJ(2));
	t22 = qJD(2) * t16 * pkin(4);
	t15 = sin(qJ(3));
	t21 = qJD(3) * t15 * pkin(5);
	t14 = -t21 - t22;
	t1 = [-t17 * t14 + (-r_i_i_C(3) * t17 + t27 * t20) * qJD(1), (t16 * t26 - t20 * t24) * pkin(4), (t15 * t26 - t20 * t23) * pkin(5), 0; t20 * t14 + (r_i_i_C(3) * t20 + t27 * t17) * qJD(1), (-t16 * t25 - t17 * t24) * pkin(4), (-t15 * t25 - t17 * t23) * pkin(5), 0; 0, -t22, -t21, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-03 01:06:52
	% EndTime: 2020-05-03 01:06:52
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (48->20), mult. (118->34), div. (0->0), fcn. (94->8), ass. (0->22)
	t89 = qJD(1) + qJD(4);
	t71 = sin(qJ(4));
	t74 = sin(qJ(1));
	t75 = cos(qJ(4));
	t78 = cos(qJ(1));
	t67 = t89 * (t71 * t74 - t75 * t78);
	t68 = t89 * (-t71 * t78 - t74 * t75);
	t88 = t67 * r_i_i_C(1) - t68 * r_i_i_C(2);
	t87 = t68 * r_i_i_C(1) + t67 * r_i_i_C(2);
	t86 = qJD(1) * t74;
	t85 = qJD(1) * t78;
	t77 = cos(qJ(2));
	t84 = qJD(2) * t77;
	t76 = cos(qJ(3));
	t83 = qJD(3) * t76;
	t73 = sin(qJ(2));
	t82 = qJD(2) * t73 * pkin(4);
	t72 = sin(qJ(3));
	t81 = qJD(3) * t72 * pkin(5);
	t70 = -t81 - t82;
	t69 = t77 * pkin(4) + t76 * pkin(5) + pkin(1) + pkin(2) + pkin(3);
	t1 = [-t69 * t85 - t70 * t74 + t88, (t73 * t86 - t78 * t84) * pkin(4), (t72 * t86 - t78 * t83) * pkin(5), t88; -t69 * t86 + t70 * t78 + t87, (-t73 * t85 - t74 * t84) * pkin(4), (-t72 * t85 - t74 * t83) * pkin(5), t87; 0, -t82, -t81, 0;];
	JaD_transl = t1;
end