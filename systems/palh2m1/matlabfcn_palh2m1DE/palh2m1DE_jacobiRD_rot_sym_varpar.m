% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-02 23:52
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh2m1DE_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh2m1DE_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh2m1DE_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_jacobiRD_rot_sym_varpar: pkin has to be [6x1] (double)');
JRD_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->7), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(2));
	t37 = qJD(2) * t31;
	t33 = cos(qJ(2));
	t36 = qJD(2) * t33;
	t35 = qJD(2) * t34;
	t30 = t32 * t37 - t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t30, t27, 0, 0; -t28, -t29, 0, 0; 0, t37, 0, 0; t29, t28, 0, 0; t27, t30, 0, 0; 0, t36, 0, 0; t39, 0, 0, 0; -t38, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t72 = qJD(2) + qJD(3);
	t74 = sin(qJ(1));
	t79 = t72 * t74;
	t75 = cos(qJ(1));
	t78 = t72 * t75;
	t77 = qJD(1) * t74;
	t76 = qJD(1) * t75;
	t73 = qJ(2) + qJ(3);
	t71 = cos(t73);
	t70 = sin(t73);
	t69 = t72 * t71;
	t68 = t72 * t70;
	t67 = t70 * t79 - t71 * t76;
	t66 = t70 * t76 + t71 * t79;
	t65 = t70 * t78 + t71 * t77;
	t64 = t70 * t77 - t71 * t78;
	t1 = [t67, t64, t64, 0; -t65, -t66, -t66, 0; 0, t68, t68, 0; t66, t65, t65, 0; t64, t67, t67, 0; 0, t69, t69, 0; t77, 0, 0, 0; -t76, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t37 = qJD(1) * sin(qJ(1));
	t36 = qJD(1) * cos(qJ(1));
	t1 = [-t36, 0, 0, 0; -t37, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; t37, 0, 0, 0; -t36, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-02 23:52:41
	% EndTime: 2020-05-02 23:52:41
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (24->5), mult. (64->6), div. (0->0), fcn. (64->4), ass. (0->8)
	t115 = qJD(1) + qJD(4);
	t110 = cos(qJ(1));
	t109 = cos(qJ(4));
	t108 = sin(qJ(1));
	t107 = sin(qJ(4));
	t102 = t115 * (-t107 * t110 - t108 * t109);
	t101 = t115 * (t107 * t108 - t109 * t110);
	t1 = [t101, 0, 0, t101; t102, 0, 0, t102; 0, 0, 0, 0; -t102, 0, 0, -t102; t101, 0, 0, t101; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
end