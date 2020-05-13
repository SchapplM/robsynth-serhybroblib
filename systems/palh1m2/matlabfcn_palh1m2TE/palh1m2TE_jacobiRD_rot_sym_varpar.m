% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% palh1m2TE
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
% pkin [22x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AM,BC,BE,BG,BL,DC,EP,GH,GP,HW,ML,OT2,T1D,T2A,T2T1,phi1,phi2,phi312,phi413,phi710,phi711]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-01 20:48
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = palh1m2TE_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh1m2TE_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh1m2TE_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'palh1m2TE_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [22 1]), ...
  'palh1m2TE_jacobiRD_rot_sym_varpar: pkin has to be [22x1] (double)');
JRD_rot=NaN(9,4);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:38
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0; -t31, 0, 0, 0; 0, 0, 0, 0; t31, 0, 0, 0; -t30, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t32 = sin(qJ(1));
	t39 = qJD(1) * t32;
	t34 = cos(qJ(1));
	t38 = qJD(1) * t34;
	t31 = sin(qJ(2));
	t37 = qJD(2) * t31;
	t33 = cos(qJ(2));
	t36 = qJD(2) * t33;
	t35 = qJD(2) * t34;
	t30 = -t32 * t37 + t33 * t38;
	t29 = t31 * t38 + t32 * t36;
	t28 = t31 * t35 + t33 * t39;
	t27 = t31 * t39 - t33 * t35;
	t1 = [t29, t28, 0, 0; t27, -t30, 0, 0; 0, -t36, 0, 0; t30, -t27, 0, 0; t28, t29, 0, 0; 0, t37, 0, 0; -t39, 0, 0, 0; t38, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->15), mult. (190->16), div. (0->0), fcn. (202->6), ass. (0->18)
	t81 = sin(qJ(3));
	t82 = sin(qJ(2));
	t84 = cos(qJ(3));
	t85 = cos(qJ(2));
	t77 = -t82 * t81 + t84 * t85;
	t92 = qJD(2) + qJD(3);
	t83 = sin(qJ(1));
	t89 = qJD(1) * t83;
	t86 = cos(qJ(1));
	t88 = qJD(1) * t86;
	t78 = t85 * t81 + t84 * t82;
	t75 = t77 * t92;
	t87 = t83 * t75 + t78 * t88;
	t72 = -t86 * t75 + t78 * t89;
	t76 = t92 * t78;
	t73 = t83 * t76 - t77 * t88;
	t71 = t86 * t76 + t77 * t89;
	t1 = [t73, t72, t72, 0; -t71, -t87, -t87, 0; 0, -t76, -t76, 0; t87, t71, t71, 0; t72, t73, t73, 0; 0, -t75, -t75, 0; -t89, 0, 0, 0; t88, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (21->6), mult. (34->14), div. (0->0), fcn. (46->8), ass. (0->14)
	t90 = qJD(1) * sin(qJ(1));
	t89 = qJD(1) * cos(qJ(1));
	t88 = cos(pkin(18));
	t86 = sin(pkin(18));
	t84 = cos(pkin(20));
	t83 = sin(pkin(20));
	t82 = pkin(22) + pkin(21);
	t81 = cos(t82);
	t80 = sin(t82);
	t79 = t86 * t83 + t88 * t84;
	t78 = -t88 * t83 + t86 * t84;
	t77 = -t78 * t81 + t80 * t79;
	t76 = t78 * t80 + t81 * t79;
	t1 = [t76 * t89, 0, 0, 0; t76 * t90, 0, 0, 0; 0, 0, 0, 0; t77 * t89, 0, 0, 0; t77 * t90, 0, 0, 0; 0, 0, 0, 0; -t90, 0, 0, 0; t89, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:40
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (127->18), mult. (256->31), div. (0->0), fcn. (300->10), ass. (0->29)
	t264 = sin(qJ(4));
	t265 = sin(qJ(1));
	t279 = t264 * t265;
	t267 = cos(qJ(4));
	t278 = t265 * t267;
	t262 = sin(pkin(20));
	t263 = cos(pkin(20));
	t266 = sin(pkin(18));
	t269 = cos(pkin(18));
	t257 = -t269 * t262 + t266 * t263;
	t258 = t266 * t262 + t269 * t263;
	t261 = pkin(22) + pkin(21);
	t259 = sin(t261);
	t260 = cos(t261);
	t256 = -t257 * t260 + t259 * t258;
	t277 = qJD(1) * t256;
	t276 = qJD(4) * t256;
	t275 = t257 * t259 + t258 * t260;
	t268 = cos(qJ(1));
	t274 = t275 * t268;
	t273 = t267 * t274 - t279;
	t272 = t264 * t268 + t275 * t278;
	t271 = t264 * t274 + t278;
	t270 = t267 * t268 - t275 * t279;
	t255 = t273 * qJD(1) + t270 * qJD(4);
	t254 = t271 * qJD(1) + t272 * qJD(4);
	t253 = t272 * qJD(1) + t271 * qJD(4);
	t252 = t270 * qJD(1) + t273 * qJD(4);
	t1 = [t255, 0, 0, t252; t253, 0, 0, t254; 0, 0, 0, -t267 * t276; -t254, 0, 0, -t253; t252, 0, 0, t255; 0, 0, 0, t264 * t276; -t268 * t277, 0, 0, 0; -t265 * t277, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (63->13), mult. (154->20), div. (0->0), fcn. (198->8), ass. (0->21)
	t84 = sin(qJ(1));
	t85 = sin(pkin(18));
	t88 = cos(pkin(17));
	t93 = sin(pkin(17));
	t94 = cos(pkin(18));
	t80 = t85 * t88 - t94 * t93;
	t81 = t93 * t85 + t88 * t94;
	t83 = sin(qJ(2));
	t86 = cos(qJ(2));
	t78 = t80 * t86 - t81 * t83;
	t89 = t78 * qJD(2);
	t90 = t83 * t80 + t81 * t86;
	t87 = cos(qJ(1));
	t91 = qJD(1) * t87;
	t101 = t84 * t89 + t90 * t91;
	t92 = qJD(1) * t84;
	t100 = -t87 * t89 + t90 * t92;
	t95 = t90 * qJD(2);
	t99 = -t78 * t91 + t84 * t95;
	t98 = t78 * t92 + t87 * t95;
	t1 = [t99, t100, 0, 0; -t98, -t101, 0, 0; 0, -t95, 0, 0; t101, t98, 0, 0; t100, t99, 0, 0; 0, -t89, 0, 0; -t92, 0, 0, 0; t91, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 7
	%% Symbolic Calculation
	% From jacobiRD_rot_7_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->3), mult. (18->10), div. (0->0), fcn. (22->6), ass. (0->9)
	t47 = qJD(1) * sin(qJ(1));
	t46 = qJD(1) * cos(qJ(1));
	t45 = cos(pkin(18));
	t43 = sin(pkin(18));
	t41 = cos(pkin(22));
	t40 = sin(pkin(22));
	t39 = t45 * t40 - t43 * t41;
	t38 = t43 * t40 + t41 * t45;
	t1 = [t38 * t46, 0, 0, 0; t38 * t47, 0, 0, 0; 0, 0, 0, 0; t39 * t46, 0, 0, 0; t39 * t47, 0, 0, 0; 0, 0, 0, 0; -t47, 0, 0, 0; t46, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 8
	%% Symbolic Calculation
	% From jacobiRD_rot_8_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:40
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (161->27), mult. (398->46), div. (0->0), fcn. (466->8), ass. (0->29)
	t167 = sin(pkin(19));
	t168 = cos(pkin(19));
	t170 = sin(qJ(2));
	t173 = cos(qJ(2));
	t163 = t167 * t173 + t168 * t170;
	t169 = sin(qJ(3));
	t172 = cos(qJ(3));
	t197 = -t167 * t170 + t168 * t173;
	t147 = t169 * t163 - t172 * t197;
	t174 = cos(qJ(1));
	t156 = t163 * qJD(2);
	t157 = t197 * qJD(2);
	t182 = t163 * t172 + t169 * t197;
	t176 = -qJD(3) * t182 - t156 * t172 - t169 * t157;
	t171 = sin(qJ(1));
	t184 = qJD(1) * t171;
	t201 = t147 * t184 + t174 * t176;
	t160 = t172 * t167 + t169 * t168;
	t198 = -t169 * t167 + t172 * t168;
	t199 = t170 * t160 - t173 * t198;
	t177 = t147 * qJD(3) + t156 * t169 - t157 * t172;
	t183 = qJD(1) * t174;
	t140 = -t171 * t177 + t182 * t183;
	t195 = t174 * t177 + t182 * t184;
	t181 = t160 * t173 + t170 * t198;
	t158 = t160 * qJD(3);
	t159 = t198 * qJD(3);
	t175 = t199 * t183 + t171 * (t181 * qJD(2) + t173 * t158 + t170 * t159);
	t1 = [t147 * t183 - t171 * t176, t195, t195, 0; t201, -t140, -t140, 0; 0, t176, t176, 0; t140, -t201, -t201, 0; t181 * t184 + t174 * (t199 * qJD(2) + t170 * t158 - t159 * t173), t175, t175, 0; 0, t177, t177, 0; -t184, 0, 0, 0; t183, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 9
	%% Symbolic Calculation
	% From jacobiRD_rot_9_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:40
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->8), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t83 = sin(qJ(1));
	t90 = qJD(1) * t83;
	t85 = cos(qJ(1));
	t89 = qJD(1) * t85;
	t82 = sin(qJ(2));
	t88 = qJD(2) * t82;
	t84 = cos(qJ(2));
	t87 = qJD(2) * t84;
	t86 = qJD(2) * t85;
	t81 = -t83 * t88 + t84 * t89;
	t80 = t82 * t89 + t83 * t87;
	t79 = t82 * t86 + t84 * t90;
	t78 = t82 * t90 - t84 * t86;
	t1 = [t80, t79, 0, 0; t78, -t81, 0, 0; 0, -t87, 0, 0; t81, -t78, 0, 0; t79, t80, 0, 0; 0, t88, 0, 0; -t90, 0, 0, 0; t89, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 10
	%% Symbolic Calculation
	% From jacobiRD_rot_10_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-05-01 20:48:39
	% EndTime: 2020-05-01 20:48:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t78 = qJ(3) + qJ(2);
	t75 = sin(t78);
	t77 = qJD(2) + qJD(3);
	t86 = t77 * t75;
	t76 = cos(t78);
	t85 = t77 * t76;
	t79 = sin(qJ(1));
	t84 = t77 * t79;
	t80 = cos(qJ(1));
	t83 = t77 * t80;
	t82 = qJD(1) * t79;
	t81 = qJD(1) * t80;
	t74 = t75 * t84 - t76 * t81;
	t73 = t75 * t81 + t76 * t84;
	t72 = t75 * t83 + t76 * t82;
	t71 = t75 * t82 - t76 * t83;
	t1 = [t74, t71, t71, 0; -t72, -t73, -t73, 0; 0, -t86, -t86, 0; t73, t72, t72, 0; t71, t74, t74, 0; 0, -t85, -t85, 0; -t82, 0, 0, 0; t81, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
end