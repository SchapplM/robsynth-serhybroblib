% Jacobian time derivative of explicit kinematic constraints of
% palh4m1IC
% with respect to passive joint coordinates
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [8x1]
%   Generalized joint coordinates (joint angles)
% qJD [8x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,CB,CE,EP,OT,TA,TD]';
% 
% Output:
% PhiD_p [(no of constraints)x(no. of passive joints)] 

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-11 23:05
% Revision: 6ae2d958c5b90587a0d08029b131cb7b66342a68 (2020-04-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function PhiD_p = palh4m1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(8,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [8 1]), ...
  'palh4m1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJ has to be [8x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [8 1]), ...
  'palh4m1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: qJD has to be [8x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'palh4m1IC_kinconstr_impl_pas_jacobianD_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_passive_jacobianD_matlab.m
% OptimizationMode: 1
% StartTime: 2020-04-11 23:05:40
% EndTime: 2020-04-11 23:05:40
% DurationCPUTime: 0.01s
% Computational Cost: add. (12->6), mult. (14->10), div. (0->0), fcn. (8->4), ass. (0->17)
unknown=NaN(3,3);
t1 = qJD(2) + qJD(4);
t2 = qJ(2) + qJ(4);
t3 = cos(t2);
t5 = t1 * t3 * pkin(2);
t6 = cos(qJ(2));
t9 = sin(qJ(2));
t12 = sin(t2);
t14 = t1 * t12 * pkin(2);
unknown(1,1) = -qJD(2) * t6 * qJ(3) - t9 * qJD(3) + t5;
unknown(1,2) = t5;
unknown(1,3) = 0.0e0;
unknown(2,1) = -qJD(2) * t9 * qJ(3) + t6 * qJD(3) + t14;
unknown(2,2) = t14;
unknown(2,3) = 0.0e0;
unknown(3,1) = 0.0e0;
unknown(3,2) = 0.0e0;
unknown(3,3) = 0.0e0;
PhiD_p  = unknown;
