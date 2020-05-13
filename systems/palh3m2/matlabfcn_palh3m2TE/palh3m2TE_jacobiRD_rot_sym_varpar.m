% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh3m2TE
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
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 01:49
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh3m2TE_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2TE_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2TE_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh3m2TE_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2TE_jacobiRD_rot_sym_varpar: pkin has to be [18x1] (double)');
JRD_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:40
	% EndTime: 2020-05-07 01:41:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:41
	% EndTime: 2020-05-07 01:41:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
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
	t1 = [t30, t27, 0, 0; -t28, -t29, 0, 0; 0, -t37, 0, 0; t29, t28, 0, 0; t27, t30, 0, 0; 0, -t36, 0, 0; -t39, 0, 0, 0; t38, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:42
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (73->12), mult. (190->16), div. (0->0), fcn. (202->6), ass. (0->18)
	t98 = qJD(2) + qJD(3);
	t92 = sin(qJ(1));
	t97 = qJD(1) * t92;
	t95 = cos(qJ(1));
	t96 = qJD(1) * t95;
	t90 = sin(qJ(3));
	t91 = sin(qJ(2));
	t93 = cos(qJ(3));
	t94 = cos(qJ(2));
	t88 = t90 * t94 + t93 * t91;
	t87 = t90 * t91 - t93 * t94;
	t85 = t98 * t88;
	t83 = -t92 * t85 - t87 * t96;
	t81 = -t95 * t85 + t87 * t97;
	t86 = t98 * t87;
	t84 = -t92 * t86 + t88 * t96;
	t82 = -t95 * t86 - t88 * t97;
	t1 = [t83, t82, t82, 0; -t81, t84, t84, 0; 0, t85, t85, 0; -t84, t81, t81, 0; t82, t83, t83, 0; 0, -t86, -t86, 0; -t97, 0, 0, 0; t96, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:42
	% EndTime: 2020-05-07 01:41:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (23->8), mult. (34->14), div. (0->0), fcn. (46->8), ass. (0->14)
	t96 = qJD(1) * sin(qJ(1));
	t95 = qJD(1) * cos(qJ(1));
	t94 = cos(pkin(15));
	t92 = sin(pkin(15));
	t90 = cos(pkin(16));
	t89 = sin(pkin(16));
	t88 = pkin(17) + pkin(18);
	t87 = cos(t88);
	t86 = sin(t88);
	t85 = -t89 * t92 + t90 * t94;
	t84 = t89 * t94 + t90 * t92;
	t83 = t86 * t84 - t85 * t87;
	t82 = t84 * t87 + t85 * t86;
	t1 = [-t83 * t95, 0, 0, 0; -t83 * t96, 0, 0, 0; 0, 0, 0, 0; t82 * t95, 0, 0, 0; t82 * t96, 0, 0, 0; 0, 0, 0, 0; -t96, 0, 0, 0; t95, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:44
	% EndTime: 2020-05-07 01:41:44
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (127->18), mult. (256->31), div. (0->0), fcn. (300->10), ass. (0->29)
	t276 = sin(qJ(4));
	t277 = sin(qJ(1));
	t291 = t276 * t277;
	t279 = cos(qJ(4));
	t290 = t277 * t279;
	t274 = sin(pkin(16));
	t275 = cos(pkin(16));
	t278 = sin(pkin(15));
	t281 = cos(pkin(15));
	t269 = t274 * t281 + t275 * t278;
	t270 = -t274 * t278 + t275 * t281;
	t273 = pkin(17) + pkin(18);
	t271 = sin(t273);
	t272 = cos(t273);
	t268 = t269 * t272 + t270 * t271;
	t289 = qJD(1) * t268;
	t288 = qJD(4) * t268;
	t287 = t269 * t271 - t270 * t272;
	t280 = cos(qJ(1));
	t286 = t287 * t280;
	t285 = -t279 * t286 - t291;
	t284 = t276 * t280 - t287 * t290;
	t283 = -t276 * t286 + t290;
	t282 = t279 * t280 + t287 * t291;
	t267 = t285 * qJD(1) + t282 * qJD(4);
	t266 = t283 * qJD(1) + t284 * qJD(4);
	t265 = t284 * qJD(1) + t283 * qJD(4);
	t264 = t282 * qJD(1) + t285 * qJD(4);
	t1 = [t267, 0, 0, t264; t265, 0, 0, t266; 0, 0, 0, -t279 * t288; -t266, 0, 0, -t265; t264, 0, 0, t267; 0, 0, 0, t276 * t288; -t280 * t289, 0, 0, 0; -t277 * t289, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:44
	% EndTime: 2020-05-07 01:41:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (63->12), mult. (154->20), div. (0->0), fcn. (198->8), ass. (0->21)
	t84 = sin(qJ(1));
	t85 = sin(pkin(15));
	t88 = cos(pkin(14));
	t93 = sin(pkin(14));
	t94 = cos(pkin(15));
	t80 = t88 * t85 - t94 * t93;
	t81 = t93 * t85 + t88 * t94;
	t83 = sin(qJ(2));
	t86 = cos(qJ(2));
	t79 = t83 * t80 - t81 * t86;
	t89 = t79 * qJD(2);
	t90 = t86 * t80 + t83 * t81;
	t87 = cos(qJ(1));
	t91 = qJD(1) * t87;
	t101 = -t84 * t89 + t90 * t91;
	t92 = qJD(1) * t84;
	t100 = t87 * t89 + t90 * t92;
	t95 = t90 * qJD(2);
	t99 = t79 * t91 + t84 * t95;
	t98 = t79 * t92 - t87 * t95;
	t1 = [t99, t100, 0, 0; t98, -t101, 0, 0; 0, -t95, 0, 0; t101, -t98, 0, 0; t100, t99, 0, 0; 0, t89, 0, 0; -t92, 0, 0, 0; t91, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiRD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->3), mult. (18->10), div. (0->0), fcn. (22->6), ass. (0->9)
	t46 = qJD(1) * sin(qJ(1));
	t45 = qJD(1) * cos(qJ(1));
	t44 = cos(pkin(15));
	t42 = sin(pkin(15));
	t40 = cos(pkin(18));
	t39 = sin(pkin(18));
	t38 = -t39 * t42 + t40 * t44;
	t37 = t39 * t44 + t40 * t42;
	t1 = [t38 * t45, 0, 0, 0; t38 * t46, 0, 0, 0; 0, 0, 0, 0; t37 * t45, 0, 0, 0; t37 * t46, 0, 0, 0; 0, 0, 0, 0; -t46, 0, 0, 0; t45, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiRD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-07 01:41:45
	% EndTime: 2020-05-07 01:41:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (57->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t87 = qJD(2) + qJD(3);
	t89 = sin(qJ(1));
	t94 = t87 * t89;
	t90 = cos(qJ(1));
	t93 = t87 * t90;
	t92 = qJD(1) * t89;
	t91 = qJD(1) * t90;
	t88 = qJ(3) + qJ(2);
	t86 = cos(t88);
	t85 = sin(t88);
	t84 = t87 * t86;
	t83 = t87 * t85;
	t82 = -t85 * t94 + t86 * t91;
	t81 = t85 * t91 + t86 * t94;
	t80 = t85 * t93 + t86 * t92;
	t79 = -t85 * t92 + t86 * t93;
	t1 = [t82, t79, t79, 0; t80, t81, t81, 0; 0, t83, t83, 0; -t81, -t80, -t80, 0; t79, t82, t82, 0; 0, t84, t84, 0; -t92, 0, 0, 0; t91, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
end