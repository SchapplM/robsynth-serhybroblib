% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% picker2Dm2OL
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [12x1]
%   Generalized joint coordinates (joint angles)
% qJD [12x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[L1,L2,L3,L4,L5,L6,e,phi05]';
% 
% Output:
% JaD_transl [3x12]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-09 23:20
% Revision: 73758128893bc0a8beabae04bd7e71472107ac81 (2020-05-07)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = picker2Dm2OL_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(12,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [12 1]), ...
  'picker2Dm2OL_jacobiaD_transl_sym_varpar: qJ has to be [12x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [12 1]), ...
  'picker2Dm2OL_jacobiaD_transl_sym_varpar: qJD has to be [12x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'picker2Dm2OL_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'picker2Dm2OL_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'picker2Dm2OL_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
JaD_transl=NaN(3,12);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(r_i_i_C(1) * t27 - r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; (r_i_i_C(1) * t26 + r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = qJD(1) + qJD(2);
	t44 = qJ(1) + qJ(2);
	t49 = sin(t44) * t43;
	t48 = cos(t44) * t43;
	t47 = r_i_i_C(1) * t49 + r_i_i_C(2) * t48;
	t46 = pkin(1) * qJD(1);
	t45 = r_i_i_C(1) * t48 - r_i_i_C(2) * t49;
	t1 = [cos(qJ(1)) * t46 + t45, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; sin(qJ(1)) * t46 + t47, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t48 = qJD(1) + qJD(2);
	t55 = pkin(2) * t48;
	t54 = pkin(1) * qJD(1);
	t49 = qJ(1) + qJ(2);
	t47 = qJ(3) + t49;
	t44 = sin(t47);
	t45 = cos(t47);
	t46 = qJD(3) + t48;
	t53 = (-r_i_i_C(1) * t45 + r_i_i_C(2) * t44) * t46;
	t52 = cos(t49) * t55 + t53;
	t51 = (-r_i_i_C(1) * t44 - r_i_i_C(2) * t45) * t46;
	t50 = sin(t49) * t55 + t51;
	t1 = [cos(qJ(1)) * t54 + t52, t52, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0; sin(qJ(1)) * t54 + t50, t50, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t52 = qJD(1) + qJD(2);
	t61 = pkin(3) * t52;
	t50 = qJD(4) + t52;
	t53 = qJ(1) + qJ(2);
	t51 = qJ(4) + t53;
	t60 = sin(t51) * t50;
	t59 = cos(t51) * t50;
	t58 = r_i_i_C(1) * t60 + r_i_i_C(2) * t59;
	t57 = pkin(1) * qJD(1);
	t56 = sin(t53) * t61 + t58;
	t55 = r_i_i_C(1) * t59 - r_i_i_C(2) * t60;
	t54 = cos(t53) * t61 + t55;
	t1 = [cos(qJ(1)) * t57 + t54, t54, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0; sin(qJ(1)) * t57 + t56, t56, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->3), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->4)
	t32 = pkin(8) + qJ(5);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(5), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, (-r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(5), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (56->6), mult. (28->9), div. (0->0), fcn. (14->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t42 = qJ(1) + qJ(2) + qJ(6);
	t39 = sin(t42);
	t40 = cos(t42);
	t41 = qJD(1) + qJD(2) + qJD(6);
	t37 = (-r_i_i_C(1) * t40 + r_i_i_C(2) * t39) * t41;
	t36 = (-r_i_i_C(1) * t39 - r_i_i_C(2) * t40) * t41;
	t1 = [cos(qJ(1)) * t43 + t37, t37, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0; sin(qJ(1)) * t43 + t36, t36, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiaD_transl_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(7));
	t26 = sin(qJ(7));
	t1 = [0, 0, 0, 0, 0, 0, (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(7), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, (r_i_i_C(1) * t27 - r_i_i_C(2) * t26) * qJD(7), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiaD_transl_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(8);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(8);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [cos(qJ(1)) * t43 + t42, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0; sin(qJ(1)) * t43 + t41, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiaD_transl_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (148->14), mult. (56->15), div. (0->0), fcn. (28->8), ass. (0->18)
	t61 = qJD(1) + qJD(2);
	t73 = pkin(2) * t61;
	t59 = qJD(3) + t61;
	t72 = pkin(6) * t59;
	t55 = qJD(9) + t59;
	t62 = qJ(1) + qJ(2);
	t60 = qJ(3) + t62;
	t58 = qJ(9) + t60;
	t71 = sin(t58) * t55;
	t70 = cos(t58) * t55;
	t69 = r_i_i_C(1) * t71 + r_i_i_C(2) * t70;
	t68 = pkin(1) * qJD(1);
	t67 = r_i_i_C(1) * t70 - r_i_i_C(2) * t71;
	t66 = -sin(t60) * t72 + t69;
	t65 = sin(t62) * t73 + t66;
	t64 = -cos(t60) * t72 + t67;
	t63 = cos(t62) * t73 + t64;
	t1 = [cos(qJ(1)) * t68 + t63, t63, t64, 0, 0, 0, 0, 0, t67, 0, 0, 0; sin(qJ(1)) * t68 + t65, t65, t66, 0, 0, 0, 0, 0, t69, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiaD_transl_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-09 23:20:39
	% EndTime: 2020-05-09 23:20:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (148->14), mult. (56->15), div. (0->0), fcn. (28->8), ass. (0->18)
	t57 = qJD(1) + qJD(2);
	t67 = pkin(3) * t57;
	t55 = qJD(4) + t57;
	t66 = pkin(4) * t55;
	t65 = pkin(1) * qJD(1);
	t58 = qJ(1) + qJ(2);
	t56 = qJ(4) + t58;
	t54 = qJ(10) + t56;
	t51 = sin(t54);
	t52 = cos(t54);
	t53 = qJD(10) + t55;
	t64 = (-r_i_i_C(1) * t52 + r_i_i_C(2) * t51) * t53;
	t63 = cos(t56) * t66 + t64;
	t62 = cos(t58) * t67 + t63;
	t61 = (-r_i_i_C(1) * t51 - r_i_i_C(2) * t52) * t53;
	t60 = sin(t56) * t66 + t61;
	t59 = sin(t58) * t67 + t60;
	t1 = [cos(qJ(1)) * t65 + t62, t62, 0, t63, 0, 0, 0, 0, 0, t64, 0, 0; sin(qJ(1)) * t65 + t59, t59, 0, t60, 0, 0, 0, 0, 0, t61, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
end